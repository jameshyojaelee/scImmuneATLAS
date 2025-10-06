"""TCR repertoire analysis using Scirpy."""

from __future__ import annotations

import json
import logging
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata as ad
import pandas as pd

from .receptor.config import (
    dataset_has_receptor,
    dataset_receptor_config,
    figures_dir as config_figures_dir,
    global_receptor_config,
    metrics_dir as config_metrics_dir,
)
from .receptor.ingest import resolved_output_path as receptor_table_path
from .receptors_io import load_receptor_table
from .utils import ensure_dir, load_config, set_seed, setup_logging, timer

logger = logging.getLogger(__name__)


def _import_scirpy():
    try:
        import scirpy as ir  # type: ignore
    except ImportError as exc:  # pragma: no cover - guarded import
        raise RuntimeError(
            "scirpy is required for TCR analysis. Please install scirpy within the analysis environment."
        ) from exc
    return ir


def _airr_from_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Convert harmonised receptor dataframe to AIRR-compliant schema for Scirpy."""
    airr_df = pd.DataFrame(
        {
            "cell_id": df["cell_id"],
            "locus": df.get("locus", df.get("chain", "")).astype(str),
            "junction": df.get("cdr3_nt", ""),
            "junction_aa": df.get("cdr3_aa", ""),
            "v_call": df.get("v_gene", ""),
            "d_call": df.get("d_gene", ""),
            "j_call": df.get("j_gene", ""),
            "productive": df.get("productive", pd.Series([pd.NA] * len(df))).astype(object),
            "duplicate_count": df.get("umis", pd.Series([pd.NA] * len(df))).astype(object),
            "consensus_count": df.get("reads", pd.Series([pd.NA] * len(df))).astype(object),
            "clonotype_id": df.get("clonotype_id", ""),
        }
    )
    return airr_df


def _merge_with_scirpy(adata: ad.AnnData, receptor_df: pd.DataFrame) -> ad.AnnData:
    ir = _import_scirpy()
    airr_df = _airr_from_dataframe(receptor_df)
    with tempfile.TemporaryDirectory() as tmpdir:
        airr_path = Path(tmpdir) / "receptor_airr.tsv"
        airr_df.to_csv(airr_path, sep="\t", index=False)
        ir.pp.merge_with_ir(adata, airr_path, how="left")
    return adata


def _compute_dataset_metrics(
    dataset_id: str,
    adata: ad.AnnData,
    receptor_df: pd.DataFrame,
    *,
    config: Dict,
) -> Dict:
    ir = _import_scirpy()

    # Define clonotypes & compute expansion/diversity
    ir.tl.define_clonotypes(adata)
    try:
        expansion = ir.tl.clonal_expansion(adata, inplace=False)
    except TypeError:  # pragma: no cover - API shim
        expansion = ir.tl.clonal_expansion(adata, target_col="clonet", inplace=False)

    expansion_counts = (
        expansion if isinstance(expansion, pd.Series) else pd.Series(expansion)
    ).fillna(0).astype(int)
    expansion_summary = expansion_counts.to_dict()

    diversity_df = ir.tl.alpha_diversity(adata, groupby=None, target_col="clonotype")
    if isinstance(diversity_df, pd.DataFrame):
        diversity = diversity_df.iloc[0].to_dict()
    else:  # pragma: no cover - Series fallback
        diversity = diversity_df.to_dict()

    # V/J usage using receptor dataframe for transparency
    v_usage = (
        receptor_df["v_gene"].replace("", pd.NA).dropna().value_counts().head(20).to_dict()
    )
    j_usage = (
        receptor_df["j_gene"].replace("", pd.NA).dropna().value_counts().head(20).to_dict()
    )

    clonotype_counts = (
        adata.obs.get("clonotype")
        .replace(["", None], pd.NA)
        .dropna()
        .value_counts()
        .to_dict()
    )

    metrics = {
        "dataset_id": dataset_id,
        "n_cells": int(adata.n_obs),
        "n_clonotypes": int(len(clonotype_counts)),
        "clonal_expansion": expansion_summary,
        "diversity": {
            key: float(value) for key, value in diversity.items() if pd.notna(value)
        },
        "v_gene_usage": v_usage,
        "j_gene_usage": j_usage,
        "top_clonotypes": sorted(
            (
                {"clonotype": clonotype, "cells": int(count)}
                for clonotype, count in clonotype_counts.items()
            ),
            key=lambda item: item["cells"],
            reverse=True,
        )[:10],
    }

    return metrics


def _load_dataset_artifacts(dataset_cfg: Dict, config: Dict) -> Tuple[ad.AnnData, pd.DataFrame]:
    dataset_id = dataset_cfg["id"]
    adata_path = Path("data/interim") / f"{dataset_id}.doublet_filtered.h5ad"
    if not adata_path.exists():
        raise FileNotFoundError(
            f"Doublet-filtered AnnData not found for dataset '{dataset_id}': {adata_path}"
        )
    adata = ad.read_h5ad(adata_path)

    receptor_table = receptor_table_path(dataset_cfg, config)
    if not receptor_table.exists():
        raise FileNotFoundError(
            f"Receptor table missing for dataset '{dataset_id}': {receptor_table}"
        )
    receptor_df = pd.read_parquet(receptor_table)
    return adata, receptor_df


def _write_json(path: Path, payload: Dict) -> None:
    ensure_dir(path.parent)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True))


def run_tcr_analysis(config: Dict) -> Dict:
    """Execute scirpy-based TCR analysis across datasets."""

    receptor_cfg = global_receptor_config(config)
    if not receptor_cfg.get("enabled", False):
        logger.info("TCR analysis disabled in configuration; skipping")
        return {}

    set_seed(config.get("seed", 0))
    metrics_dir = config_metrics_dir(config)
    ensure_dir(metrics_dir)

    datasets = [
        dataset_cfg
        for dataset_cfg in config.get("datasets", [])
        if dataset_has_receptor(dataset_cfg)
    ]

    if not datasets:
        logger.warning("TCR analysis requested but no datasets provide receptor data")
        return {}

    dataset_metrics: Dict[str, Dict] = {}
    aggregated_receptor: List[pd.DataFrame] = []

    for dataset_cfg in datasets:
        dataset_id = dataset_cfg["id"]
        try:
            with timer(f"TCR analysis ({dataset_id})"):
                adata, receptor_df = _load_dataset_artifacts(dataset_cfg, config)
                aggregated_receptor.append(receptor_df.assign(dataset_id=dataset_id))
                _merge_with_scirpy(adata, receptor_df)
                metrics = _compute_dataset_metrics(dataset_id, adata, receptor_df, config=config)
                dataset_metrics[dataset_id] = metrics
                _write_json(metrics_dir / f"{dataset_id}_tcr_metrics.json", metrics)
        except Exception as exc:  # pragma: no cover - runtime safety
            logger.exception("Failed TCR analysis for %s: %s", dataset_id, exc)
            raise

    aggregated_df = pd.concat(aggregated_receptor, ignore_index=True)
    global_metrics = {
        "total_cells": int(sum(m["n_cells"] for m in dataset_metrics.values())),
        "total_clonotypes": int(
            aggregated_df["clonotype_id"].replace("", pd.NA).dropna().nunique()
        ),
        "v_gene_usage": (
            aggregated_df["v_gene"].replace("", pd.NA).dropna().value_counts().head(25).to_dict()
        ),
        "j_gene_usage": (
            aggregated_df["j_gene"].replace("", pd.NA).dropna().value_counts().head(25).to_dict()
        ),
    }

    summary = {"datasets": dataset_metrics, "global": global_metrics}
    _write_json(metrics_dir / "tcr_summary.json", summary)
    logger.info("TCR analysis complete for %d datasets", len(dataset_metrics))
    return summary


def main() -> None:  # pragma: no cover - CLI helper
    import argparse

    parser = argparse.ArgumentParser(description="Run TCR repertoire analysis with Scirpy")
    parser.add_argument("--config", default="config/atlas.yaml", help="YAML configuration path")
    parser.add_argument("--log-level", default="INFO", help="Logging level (default: INFO)")
    args = parser.parse_args()

    setup_logging(args.log_level)
    config = load_config(args.config)
    run_tcr_analysis(config)


if __name__ == "__main__":  # pragma: no cover
    main()
