"""Cell type annotation functions for the Single-cell Immune Atlas."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from .utils import ensure_dir, load_config, setup_logging, timer


def load_immune_markers() -> pd.DataFrame:
    """Load immune cell type markers."""
    markers_path = Path(__file__).parent / "markers" / "immune_markers_human.tsv"

    if not markers_path.exists():  # pragma: no cover - depends on repository state
        raise FileNotFoundError(f"Markers file not found: {markers_path}")

    markers_df = pd.read_csv(markers_path, sep="\t")

    required_cols = ["cell_type", "gene", "direction"]
    missing = [col for col in required_cols if col not in markers_df.columns]
    if missing:
        raise ValueError(f"Markers file must contain columns: {required_cols}")

    markers_df["direction"] = markers_df["direction"].str.lower()
    logging.info(
        "Loaded %d markers covering %d cell types",
        len(markers_df),
        markers_df["cell_type"].nunique(),
    )
    return markers_df


def _augment_var_with_feature_names_from_interim(adata: ad.AnnData) -> None:
    """Populate `feature_name` for integrated data if available in interim outputs."""
    if "feature_name" in adata.var.columns:
        return

    interim_dir = Path("data/interim")
    if not interim_dir.exists():  # pragma: no cover - depends on runtime execution
        return

    mapping = {}
    for h5ad_path in sorted(interim_dir.glob("*.h5ad")):
        try:
            ds = ad.read_h5ad(h5ad_path)
        except Exception:  # pragma: no cover - defensive
            continue
        var = ds.var
        if "feature_name" not in var.columns:
            continue
        if "feature_id" in var.columns:
            mapping.update(dict(zip(var["feature_id"].astype(str), var["feature_name"].astype(str))))
        if "soma_joinid" in var.columns:
            mapping.update(dict(zip(var["soma_joinid"].astype(str), var["feature_name"].astype(str))))

    if mapping:
        adata.var["feature_name"] = [mapping.get(str(v), str(v)) for v in adata.var_names]


def _select_marker_symbols(adata: ad.AnnData) -> pd.Series:
    if "feature_name" in adata.var.columns:
        return adata.var["feature_name"].astype(str)
    if getattr(adata, "raw", None) is not None and adata.raw is not None and "feature_name" in adata.raw.var.columns:
        raw_map = adata.raw.var["feature_name"].astype(str)
        mapped = [raw_map.get(v, v) for v in adata.var_names.astype(str)]
        return pd.Series(mapped, index=adata.var_names, dtype=str)
    return pd.Series(adata.var_names.astype(str), index=adata.var_names, dtype=str)


def _score_cell_type(
    adata: ad.AnnData,
    positive_genes: Iterable[str],
    negative_genes: Iterable[str],
    score_name: str,
    ctrl_size: int = 50,
) -> np.ndarray:
    """Score a cell type using positive/negative markers with background control."""
    positive_genes = list(dict.fromkeys(g for g in positive_genes if g))
    negative_genes = list(dict.fromkeys(g for g in negative_genes if g))

    if not positive_genes and not negative_genes:
        return np.zeros(adata.n_obs)

    obs_key_pos = f"{score_name}_pos"
    obs_key_neg = f"{score_name}_neg"

    if positive_genes:
        sc.tl.score_genes(
            adata,
            gene_list=positive_genes,
            ctrl_size=min(ctrl_size, max(len(positive_genes), 1)),
            score_name=obs_key_pos,
            use_raw=False,
        )
        pos_scores = adata.obs[obs_key_pos].to_numpy()
    else:
        pos_scores = np.zeros(adata.n_obs)

    if negative_genes:
        sc.tl.score_genes(
            adata,
            gene_list=negative_genes,
            ctrl_size=min(ctrl_size, max(len(negative_genes), 1)),
            score_name=obs_key_neg,
            use_raw=False,
        )
        neg_scores = adata.obs[obs_key_neg].to_numpy()
    else:
        neg_scores = np.zeros(adata.n_obs)

    # Clean up intermediate columns
    for col in (obs_key_pos, obs_key_neg):
        if col in adata.obs:
            del adata.obs[col]

    return pos_scores - neg_scores


def compute_marker_scores(adata: ad.AnnData, markers_df: pd.DataFrame) -> Dict[str, np.ndarray]:
    """Compute marker-based scores for each cell type using gene-signature scoring."""
    symbol_series = _select_marker_symbols(adata)

    scores: Dict[str, np.ndarray] = {}
    for cell_type, subset in markers_df.groupby("cell_type"):
        pos_genes = subset.loc[subset["direction"] != "down", "gene"].astype(str).tolist()
        neg_genes = subset.loc[subset["direction"] == "down", "gene"].astype(str).tolist()

        pos_genes = [g for g in pos_genes if g in symbol_series.values]
        neg_genes = [g for g in neg_genes if g in symbol_series.values]

        if not pos_genes and not neg_genes:
            logging.warning("No markers matched for %s", cell_type)
            scores[cell_type] = np.zeros(adata.n_obs)
            continue

        scores[cell_type] = _score_cell_type(
            adata,
            positive_genes=pos_genes,
            negative_genes=neg_genes,
            score_name=f"{cell_type}_score_tmp",
        )
        logging.info(
            "%s: %d up markers, %d down markers",
            cell_type,
            len(pos_genes),
            len(neg_genes),
        )

    return scores


def _apply_confidence_thresholds(
    scores_df: pd.DataFrame,
    percentile: float,
) -> pd.Series:
    """Determine confident assignments based on percentile thresholds per cell type."""
    best_matches = scores_df.idxmax(axis=1)
    best_scores = scores_df.max(axis=1)

    thresholds = {
        cell_type: np.quantile(scores_df[cell_type].dropna(), percentile)
        if not scores_df[cell_type].dropna().empty
        else np.inf
        for cell_type in scores_df.columns
    }

    confident = best_matches.copy()
    for idx, cell_type in best_matches.items():
        score = best_scores.loc[idx]
        threshold = thresholds.get(cell_type, np.inf)
        if np.isnan(score) or score < threshold:
            confident.loc[idx] = "Unknown"

    return confident


def _write_annotation_metrics(
    adata: ad.AnnData,
    outputs_cfg: Dict,
    original_labels: Optional[pd.Series],
    scores_df: pd.DataFrame,
) -> None:
    metrics_dir = Path(outputs_cfg.get("metrics_dir", "processed/metrics"))
    ensure_dir(metrics_dir)

    # Summary of final assignments
    summary = adata.obs["cell_type"].value_counts().to_dict()
    summary_path = metrics_dir / "annotation_summary.json"
    with open(summary_path, "w") as f:
        json.dump({"cell_type_counts": summary}, f, indent=2)
    logging.info("Saved annotation summary to %s", summary_path)

    # Save per-cell scores for downstream QC
    scores_path = metrics_dir / "annotation_scores.tsv"
    scores_df.to_csv(scores_path, sep="\t")
    logging.info("Saved annotation scores to %s", scores_path)

    # Confusion matrix versus prior labels if available
    confusion_path = metrics_dir / "annotation_confusion.tsv"
    if original_labels is not None and not original_labels.isna().all():
        confusion = pd.crosstab(
            original_labels.fillna("Unlabeled"),
            adata.obs["cell_type_pred"].astype(str),
            dropna=False,
        )
        confusion.to_csv(confusion_path, sep="\t")
        logging.info("Saved annotation confusion matrix to %s", confusion_path)
    else:
        confusion_path.write_text("# No reference labels available for confusion matrix\n")


def score_and_annotate(
    adata: ad.AnnData,
    annotation_cfg: Dict,
    outputs_cfg: Dict,
) -> ad.AnnData:
    """Score markers and annotate cell types."""
    logging.info("Starting marker-based annotation")

    adata = adata.copy()
    original_labels = adata.obs["cell_type"].copy() if "cell_type" in adata.obs.columns else None
    adata.obs["cell_type_initial"] = original_labels

    if "feature_name" not in adata.var.columns:
        _augment_var_with_feature_names_from_interim(adata)

    markers_df = load_immune_markers()
    scores = compute_marker_scores(adata, markers_df)

    for cell_type, values in scores.items():
        adata.obs[f"{cell_type}_score"] = values

    scores_df = pd.DataFrame(scores, index=adata.obs_names)
    confidence_pct = annotation_cfg.get("min_confidence_pct", 0.75)
    predictions = _apply_confidence_thresholds(scores_df, confidence_pct)

    adata.obs["cell_type_pred"] = pd.Categorical(predictions.astype(str))

    # Merge with existing labels if provided
    if original_labels is None or original_labels.isna().all():
        adata.obs["cell_type"] = adata.obs["cell_type_pred"]
    else:
        updated = original_labels.astype(str).copy()
        unknown_mask = updated.isna() | (updated.str.lower() == "unknown")
        updated.loc[unknown_mask] = predictions.loc[unknown_mask]
        adata.obs["cell_type"] = pd.Categorical(updated.fillna("Unknown"))

    # Optional reference mapping fallback
    reference_cfg = annotation_cfg.get("reference_mapping") or {}
    label_key = reference_cfg.get("label_key")
    if label_key and label_key in adata.obs:
        unknown_mask = adata.obs["cell_type"].astype(str) == "Unknown"
        fallback = adata.obs[label_key].astype(str)
        adata.obs.loc[unknown_mask & fallback.notna(), "cell_type"] = fallback[unknown_mask & fallback.notna()]
        logging.info("Applied reference mapping fallback from %s", label_key)

    compartment_map = {
        "CD8_T": "T_cell",
        "CD4_T": "T_cell",
        "Treg": "T_cell",
        "NK": "NK_cell",
        "B_cell": "B_cell",
        "Plasma": "B_cell",
        "Mono": "Myeloid",
        "DC_cDC1": "Myeloid",
        "DC_cDC2": "Myeloid",
        "Unknown": "Unknown",
    }
    adata.obs["compartment"] = (
        adata.obs["cell_type"].astype(str).map(compartment_map).fillna("Unknown")
    )

    _write_annotation_metrics(adata, outputs_cfg, original_labels, scores_df)

    logging.info("Annotation completed")
    return adata


def run_annotation(config: Dict) -> None:
    """Run cell type annotation."""
    input_path = Path("processed") / "integrated_atlas.h5ad"
    logging.info("Loading integrated data from %s", input_path)
    adata = ad.read_h5ad(input_path)

    annotation_cfg = config.get("annotation", {})
    outputs_cfg = config.get("outputs", {})
    adata_annotated = score_and_annotate(adata, annotation_cfg, outputs_cfg)

    output_path = Path("processed") / "integrated_annotated.h5ad"
    adata_annotated.write(output_path)
    logging.info("Saved annotated atlas to %s", output_path)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run cell type annotation")
    parser.add_argument("--config", default="config/atlas.yaml", help="Config file path")

    args = parser.parse_args()

    setup_logging()
    cfg = load_config(args.config)

    with timer("Cell type annotation"):
        run_annotation(cfg)
