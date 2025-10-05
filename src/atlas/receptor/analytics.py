"""Repertoire-level analytics for immune receptor data."""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import anndata as ad
import numpy as np
import pandas as pd

from ..utils import ensure_dir
from . import ingest
from .config import (
    dataset_has_receptor,
    global_receptor_config,
    metrics_dir as config_metrics_dir,
    figures_dir as config_figures_dir,
)

SUMMARY_FILENAME = "repertoire_summary.json"
EXPANSION_FILENAME = "clonotype_expansion.tsv"
DIVERSITY_FILENAME = "diversity_indices.json"
VJ_USAGE_FILENAME = "vj_usage.tsv"
PAIRING_FILENAME = "chain_pairing_summary.json"


@dataclass
class AnalyticsResult:
    combined: pd.DataFrame
    cell_summary: pd.DataFrame
    clonotype_summary: pd.DataFrame
    metrics_dir: Path
    figures_dir: Path
    summary_path: Path
    expansion_path: Path
    diversity_path: Path
    vj_usage_path: Path
    pairing_path: Path
    summary_payload: Dict


def _dataset_frames(config: Dict) -> Dict[str, pd.DataFrame]:
    frames: Dict[str, pd.DataFrame] = {}
    for dataset_cfg in config.get("datasets", []):
        if not (dataset_has_receptor(dataset_cfg) or dataset_cfg.get("receptor")):
            continue
        dataset_id = dataset_cfg["id"]
        parquet_path = ingest.resolved_output_path(dataset_cfg, config)
        if not parquet_path.exists():
            logging.warning("Skipping %s receptor analytics (missing %s)", dataset_id, parquet_path)
            continue
        frames[dataset_id] = pd.read_parquet(parquet_path)
    return frames


def _filter_productive(df: pd.DataFrame, productive_only: bool) -> pd.DataFrame:
    if not productive_only:
        return df
    if "productive" not in df.columns:
        return df
    productive_mask = df["productive"].fillna(False)
    return df[productive_mask]


def _default_expansion_labels(n_bins: int) -> List[str]:
    base = ["singleton", "small", "intermediate", "expanded"]
    if n_bins <= len(base):
        return base[:n_bins]
    labels = base[:]
    for idx in range(len(base), n_bins):
        labels.append(f"bin_{idx}")
    return labels


def _assign_expansion(size: int, thresholds: List[int], labels: List[str]) -> str:
    for idx, threshold in enumerate(thresholds):
        if size <= threshold:
            return labels[idx]
    return labels[len(thresholds)]


def _clonotype_summary(df: pd.DataFrame, thresholds: List[int], labels: List[str]) -> pd.DataFrame:
    usable = df[df["clonotype_id"].replace("", pd.NA).notna()].copy()
    if usable.empty:
        return pd.DataFrame(columns=["dataset_id", "clonotype_id", "cell_count", "expansion_bin"])
    grouped = (
        usable.groupby(["dataset_id", "clonotype_id"])["cell_id"].nunique().reset_index(name="cell_count")
    )
    grouped["expansion_bin"] = grouped["cell_count"].apply(
        lambda size: _assign_expansion(int(size), thresholds, labels)
    )
    return grouped.sort_values(["dataset_id", "cell_count"], ascending=[True, False]).reset_index(drop=True)


def _paired_label(chain_loci: set[str]) -> str:
    loci = {loc for loc in chain_loci if loc}
    if {"TRA", "TRB"}.issubset(loci):
        return "TCR_alpha_beta"
    if {"TRG", "TRD"}.issubset(loci):
        return "TCR_gamma_delta"
    if "IGH" in loci and ("IGK" in loci or "IGL" in loci):
        return "BCR_heavy_light"
    return "unpaired"


def _cell_summary(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=[
            "dataset_id",
            "cell_id",
            "clonotype_id",
            "n_chains",
            "n_productive_chains",
            "has_productive_chain",
            "paired_label",
        ])
    grouped = df.groupby("cell_id")
    rows = []
    for cell_id, subset in grouped:
        loci = set(subset.get("locus", pd.Series(dtype=str)).dropna().astype(str))
        productive = subset.get("productive", pd.Series(dtype=bool)).fillna(False)
        clonotypes = subset.get("clonotype_id", pd.Series(dtype=str)).replace("", pd.NA).dropna().unique()
        rows.append({
            "dataset_id": subset["dataset_id"].iloc[0],
            "cell_id": cell_id,
            "clonotype_id": clonotypes[0] if len(clonotypes) else "",
            "n_chains": int(len(subset)),
            "n_productive_chains": int(productive.sum()),
            "has_productive_chain": bool(productive.any()),
            "paired_label": _paired_label(loci),
        })
    return pd.DataFrame(rows)


def _diversity_metrics(counts: pd.Series) -> Dict[str, float]:
    total = counts.sum()
    if total == 0:
        return {"shannon": 0.0, "simpson": 0.0, "gini": 0.0}
    proportions = counts / total
    shannon = float(-(proportions * np.log(proportions + 1e-12)).sum())
    simpson = float(1.0 - (proportions.pow(2)).sum())
    sorted_counts = np.sort(counts.values)
    n = len(sorted_counts)
    if n == 0 or total == 0:
        gini = 0.0
    else:
        index = np.arange(1, n + 1)
        gini = float((np.sum((2 * index - n - 1) * sorted_counts)) / (n * total))
    return {"shannon": shannon, "simpson": simpson, "gini": gini}


def _write_json(path: Path, payload: Dict) -> None:
    ensure_dir(path.parent)
    path.write_text(json.dumps(payload, indent=2))


def _update_atlas_obs(config: Dict, cell_summary: pd.DataFrame, clonotype_summary: pd.DataFrame) -> None:
    atlas_path = Path(global_receptor_config(config).get("atlas_path"))
    if not atlas_path.exists():
        logging.warning("Annotated atlas not found for receptor augmentation: %s", atlas_path)
        return

    adata = ad.read_h5ad(atlas_path)
    per_cell = cell_summary.set_index("cell_id") if not cell_summary.empty else pd.DataFrame()
    if not per_cell.empty:
        reindexed = per_cell.reindex(adata.obs_names)
        for column in ["clonotype_id", "n_chains", "n_productive_chains", "paired_label", "has_productive_chain"]:
            adata.obs[column] = reindexed[column]
    else:
        for column in ["clonotype_id", "n_chains", "n_productive_chains", "paired_label", "has_productive_chain"]:
            adata.obs[column] = None

    if not clonotype_summary.empty:
        size_map = clonotype_summary.set_index("clonotype_id")["cell_count"].to_dict()
        bin_map = clonotype_summary.set_index("clonotype_id")["expansion_bin"].to_dict()
        adata.obs["clonotype_size"] = adata.obs["clonotype_id"].map(size_map)
        adata.obs["clonotype_expansion_bin"] = adata.obs["clonotype_id"].map(bin_map)
    else:
        adata.obs["clonotype_size"] = None
        adata.obs["clonotype_expansion_bin"] = None

    adata.write(atlas_path)
    logging.info("Annotated atlas updated with receptor annotations (%s)", atlas_path)


def run_analytics(config: Dict) -> AnalyticsResult:
    frames = _dataset_frames(config)
    if not frames:
        logging.warning("No receptor tables available for analytics")
        metrics_dir = config_metrics_dir(config)
        figures_dir = config_figures_dir(config)
        ensure_dir(metrics_dir)
        ensure_dir(figures_dir)
        summary_path = metrics_dir / SUMMARY_FILENAME
        expansion_path = metrics_dir / EXPANSION_FILENAME
        diversity_path = metrics_dir / DIVERSITY_FILENAME
        vj_usage_path = metrics_dir / VJ_USAGE_FILENAME
        pairing_path = metrics_dir / PAIRING_FILENAME
        _write_json(summary_path, {"datasets": {}, "global": {}})
        _write_json(diversity_path, {})
        expansion_path.write_text("dataset_id\tclonotype_id\tcell_count\texpansion_bin\n")
        vj_usage_path.write_text("dataset_id\tv_gene\tj_gene\tcount\n")
        _write_json(pairing_path, {})
        return AnalyticsResult(
            combined=pd.DataFrame(),
            cell_summary=pd.DataFrame(),
            clonotype_summary=pd.DataFrame(),
            metrics_dir=metrics_dir,
            figures_dir=figures_dir,
            summary_path=summary_path,
            expansion_path=expansion_path,
            diversity_path=diversity_path,
            vj_usage_path=vj_usage_path,
            pairing_path=pairing_path,
            summary_payload={"datasets": {}, "global": {}},
        )

    receptor_cfg = global_receptor_config(config)
    productive_only = bool(receptor_cfg.get("productive_only", False))
    thresholds = sorted(receptor_cfg.get("expansion_bins", [1, 5, 20]))
    labels = receptor_cfg.get("expansion_labels")
    if not labels or len(labels) != len(thresholds) + 1:
        labels = _default_expansion_labels(len(thresholds) + 1)

    datasets_combined: List[pd.DataFrame] = []
    per_cell_frames: List[pd.DataFrame] = []
    dataset_payload: Dict[str, Dict] = {}

    for dataset_id, frame in frames.items():
        filtered = _filter_productive(frame, productive_only)
        datasets_combined.append(filtered.assign(dataset_id=dataset_id))
        cell_metrics = _cell_summary(filtered)
        per_cell_frames.append(cell_metrics)

    combined = pd.concat(datasets_combined, ignore_index=True) if datasets_combined else pd.DataFrame()
    cell_summary = pd.concat(per_cell_frames, ignore_index=True) if per_cell_frames else pd.DataFrame()

    clonotype_summary = _clonotype_summary(combined, thresholds, labels) if not combined.empty else pd.DataFrame()

    metrics_dir = config_metrics_dir(config)
    figures_dir = config_figures_dir(config)
    ensure_dir(metrics_dir)
    ensure_dir(figures_dir)

    # Diversity per dataset
    diversity_payload: Dict[str, Dict[str, float]] = {}
    for dataset_id, df_dataset in combined.groupby("dataset_id"):
        counts = df_dataset[df_dataset["clonotype_id"].replace("", pd.NA).notna()].groupby("clonotype_id")[
            "cell_id"
        ].nunique()
        diversity_payload[dataset_id] = _diversity_metrics(counts)
        top_clones = clonotype_summary[clonotype_summary["dataset_id"] == dataset_id].head(10)
        dataset_payload[dataset_id] = {
            "n_cells": int(df_dataset["cell_id"].nunique()),
            "n_sequences": int(len(df_dataset)),
            "n_clonotypes": int(counts.shape[0]),
            "diversity": diversity_payload[dataset_id],
            "top_clonotypes": top_clones.to_dict(orient="records"),
        }

    global_counts = combined[combined["clonotype_id"].replace("", pd.NA).notna()].groupby("clonotype_id")[
        "cell_id"
    ].nunique()
    global_summary = {
        "n_cells": int(combined["cell_id"].nunique()),
        "n_sequences": int(len(combined)),
        "n_clonotypes": int(global_counts.shape[0]),
        "diversity": _diversity_metrics(global_counts),
    }

    summary_payload = {"datasets": dataset_payload, "global": global_summary}

    summary_path = metrics_dir / SUMMARY_FILENAME
    expansion_path = metrics_dir / EXPANSION_FILENAME
    diversity_path = metrics_dir / DIVERSITY_FILENAME
    vj_usage_path = metrics_dir / VJ_USAGE_FILENAME
    pairing_path = metrics_dir / PAIRING_FILENAME

    _write_json(summary_path, summary_payload)
    _write_json(diversity_path, diversity_payload)
    clonotype_summary.to_csv(expansion_path, sep="\t", index=False)

    # V/J usage
    if not combined.empty:
        v_usage = (
            combined[combined["v_gene"] != ""].groupby(["dataset_id", "v_gene"]).size().reset_index(name="count")
        )
        if "j_gene" in combined.columns:
            vj_usage = (
                combined[(combined["v_gene"] != "") & (combined["j_gene"] != "")]
                .groupby(["dataset_id", "v_gene", "j_gene"])
                .size()
                .reset_index(name="count")
            )
        else:
            vj_usage = pd.DataFrame(columns=["dataset_id", "v_gene", "j_gene", "count"])
    else:
        v_usage = pd.DataFrame(columns=["dataset_id", "v_gene", "count"])
        vj_usage = pd.DataFrame(columns=["dataset_id", "v_gene", "j_gene", "count"])
    vj_usage.to_csv(vj_usage_path, sep="\t", index=False)

    pairing_counts = (
        cell_summary.groupby(["dataset_id", "paired_label"]).size().reset_index(name="count")
        if not cell_summary.empty
        else pd.DataFrame(columns=["dataset_id", "paired_label", "count"])
    )
    pairing_payload = {
        dataset_id: group.set_index("paired_label")["count"].to_dict()
        for dataset_id, group in pairing_counts.groupby("dataset_id")
    }
    _write_json(pairing_path, pairing_payload)

    _write_json(metrics_dir / "v_usage_summary.json", v_usage.to_dict(orient="records"))

    _update_atlas_obs(config, cell_summary, clonotype_summary)

    return AnalyticsResult(
        combined=combined,
        cell_summary=cell_summary,
        clonotype_summary=clonotype_summary,
        metrics_dir=metrics_dir,
        figures_dir=figures_dir,
        summary_path=summary_path,
        expansion_path=expansion_path,
        diversity_path=diversity_path,
        vj_usage_path=vj_usage_path,
        pairing_path=pairing_path,
        summary_payload=summary_payload,
    )
