"""Quality control functions for the Single-cell Immune Atlas."""

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from .io import load_matrix
from .utils import ensure_dir, load_config, setup_logging, timer


def compute_qc_metrics(adata: ad.AnnData) -> ad.AnnData:
    """Compute QC metrics including mitochondrial gene percentage."""
    # Make a copy to avoid modifying original
    adata = adata.copy()
    
    # Identify mitochondrial genes (Scanpy convention: genes prefixed with MT-)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    
    # Calculate QC metrics, including mitochondrial stats
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    
    # Ensure pct_counts_mt exists (older versions may not add it)
    if "pct_counts_mt" not in adata.obs.columns and "total_counts_mt" in adata.obs.columns:
        adata.obs["pct_counts_mt"] = (
            adata.obs["total_counts_mt"] / adata.obs["total_counts"] * 100
        )
    
    logging.info(f"QC metrics computed for {adata.n_obs} cells, {adata.n_vars} genes")
    logging.info(
        "Mean genes per cell: %0.1f",
        adata.obs.get("n_genes", adata.obs.get("n_genes_by_counts")).mean(),
    )
    logging.info(f"Mean MT%: {adata.obs['pct_counts_mt'].mean():.2f}%")

    return adata


def _qc_pass_mask(adata: ad.AnnData, qc_config: Dict) -> np.ndarray:
    """Compute boolean mask indicating cells that pass QC thresholds."""
    obs = adata.obs
    mask = np.ones(adata.n_obs, dtype=bool)

    min_genes = qc_config.get("min_genes")
    if min_genes is not None:
        metric = obs.get("n_genes_by_counts", obs.get("n_genes"))
        mask &= metric >= min_genes

    max_genes = qc_config.get("max_genes")
    if max_genes is not None:
        metric = obs.get("n_genes_by_counts", obs.get("n_genes"))
        mask &= metric <= max_genes

    min_counts = qc_config.get("min_counts")
    if min_counts is not None:
        mask &= obs["total_counts"] >= min_counts

    max_counts = qc_config.get("max_counts")
    if max_counts is not None:
        mask &= obs["total_counts"] <= max_counts

    max_mt_pct = qc_config.get("max_mt_pct")
    if max_mt_pct is not None:
        mask &= obs["pct_counts_mt"] < max_mt_pct

    return mask


def apply_filters(adata: ad.AnnData, qc_config: Dict) -> Tuple[ad.AnnData, np.ndarray, Dict[str, float]]:
    """Apply QC filters based on configuration and return filtered AnnData, mask, summary."""
    adata = adata.copy()
    mask = _qc_pass_mask(adata, qc_config)

    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    filtered = adata[mask].copy()

    # Filter genes (keep genes expressed in at least min_cells_per_gene)
    min_cells_per_gene = qc_config.get("min_cells_per_gene", 3)
    sc.pp.filter_genes(filtered, min_cells=min_cells_per_gene)

    n_cells_after = filtered.n_obs
    n_genes_after = filtered.n_vars

    logging.info(
        "QC filtering: %d → %d cells (%.1f%% retained)",
        n_cells_before,
        n_cells_after,
        100 * n_cells_after / max(n_cells_before, 1),
    )
    logging.info(
        "QC filtering: %d → %d genes (%.1f%% retained)",
        n_genes_before,
        n_genes_after,
        100 * n_genes_after / max(n_genes_before, 1),
    )

    summary = {
        "n_cells_before": int(n_cells_before),
        "n_cells_after": int(n_cells_after),
        "n_genes_before": int(n_genes_before),
        "n_genes_after": int(n_genes_after),
        "pct_cells_retained": float(100 * n_cells_after / max(n_cells_before, 1)),
        "pct_genes_retained": float(100 * n_genes_after / max(n_genes_before, 1)),
    }

    return filtered, mask, summary


def _write_qc_tables(
    dataset_id: str,
    adata: ad.AnnData,
    mask: np.ndarray,
    qc_config: Dict,
    summary: Dict[str, float],
) -> None:
    """Persist per-cell QC metrics and dataset summary to disk."""
    metrics_dir = Path(qc_config.get("metrics_dir", "processed/metrics"))
    ensure_dir(metrics_dir)

    qc_table = adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy()
    qc_table["qc_pass"] = mask
    table_path = Path("data/interim") / f"{dataset_id}_qc.tsv"
    qc_table.to_csv(table_path, sep="\t")
    logging.info("Saved QC table to %s", table_path)

    summary_path = metrics_dir / f"{dataset_id}_qc_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    logging.info("Saved QC summary to %s", summary_path)


def _write_qc_plots(adata: ad.AnnData, dataset_id: str, qc_config: Dict) -> None:
    """Generate violin plots for QC metrics."""
    plots_dir = Path(qc_config.get("qc_plots_dir", "processed/figures/qc"))
    ensure_dir(plots_dir)

    sample_df = adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy()
    sample_size = min(5000, sample_df.shape[0])
    if sample_df.shape[0] > sample_size:
        sample_df = sample_df.sample(sample_size, random_state=0)

    if sample_df.empty:
        for suffix in ("umi_violin", "mt_violin"):
            fig, ax = plt.subplots(figsize=(4, 3))
            ax.axis("off")
            ax.text(0.5, 0.5, "No data available", ha="center", va="center")
            placeholder = plots_dir / f"{dataset_id}_{suffix}.png"
            fig.savefig(placeholder, dpi=300)
            plt.close(fig)
        logging.warning("No cells available for QC plotting for %s", dataset_id)
        return

    # UMI / gene counts violin plots
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    log_counts = np.log10(sample_df["total_counts"] + 1)
    sns.violinplot(y=log_counts, ax=axes[0], color="#4c72b0")
    axes[0].set_title("log10 Total Counts")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("log10(counts+1)")

    sns.violinplot(y=sample_df["n_genes_by_counts"], ax=axes[1], color="#55a868")
    axes[1].set_title("Detected Genes")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("Genes")

    fig.suptitle(f"QC distributions – {dataset_id}")
    plt.tight_layout()
    umi_path = plots_dir / f"{dataset_id}_umi_violin.png"
    fig.savefig(umi_path, dpi=300)
    plt.close(fig)
    logging.info("Saved QC violin plot to %s", umi_path)

    # Mitochondrial percentage violin plot
    fig, ax = plt.subplots(figsize=(5, 4))
    sns.violinplot(y=sample_df["pct_counts_mt"], ax=ax, color="#c44e52")
    ax.set_title("Mitochondrial %")
    ax.set_xlabel("")
    ax.set_ylabel("pct_counts_mt")
    fig.suptitle(f"QC distributions – {dataset_id}")
    plt.tight_layout()
    mt_path = plots_dir / f"{dataset_id}_mt_violin.png"
    fig.savefig(mt_path, dpi=300)
    plt.close(fig)
    logging.info("Saved QC MT plot to %s", mt_path)


def qc_summary_table(adata_list: List[ad.AnnData]) -> pd.DataFrame:
    """Generate a summary table of QC metrics across datasets.

    Falls back to NaN when optional QC columns are unavailable and returns an empty
    DataFrame when no inputs are provided.
    """

    def _dataset_label(adata: ad.AnnData) -> str:
        obs = adata.obs
        if "dataset_id" in obs.columns and not obs["dataset_id"].empty:
            return str(obs["dataset_id"].iloc[0])
        return str(adata.uns.get("dataset_id", "unknown"))

    def _obs_stat(adata: ad.AnnData, column: str, reducer: str) -> float:
        obs = adata.obs
        if column not in obs.columns or obs.empty:
            return float("nan")
        series = obs[column]
        if series.empty:
            return float("nan")
        func = getattr(series, reducer)
        try:
            return float(func())
        except Exception:  # pragma: no cover - defensive
            return float("nan")

    rows: List[Dict[str, Any]] = []
    for adata in adata_list:
        if not isinstance(adata, ad.AnnData):
            continue

        row = {
            "dataset_id": _dataset_label(adata),
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "mean_genes_per_cell": _obs_stat(adata, "n_genes_by_counts", "mean"),
            "median_genes_per_cell": _obs_stat(adata, "n_genes_by_counts", "median"),
            "mean_counts_per_cell": _obs_stat(adata, "total_counts", "mean"),
            "median_counts_per_cell": _obs_stat(adata, "total_counts", "median"),
            "mean_mt_pct": _obs_stat(adata, "pct_counts_mt", "mean"),
            "median_mt_pct": _obs_stat(adata, "pct_counts_mt", "median"),
        }
        rows.append(row)

    return pd.DataFrame(rows)


def process_dataset_qc(dataset_id: str, config: Dict) -> None:
    """Process QC for a single dataset."""
    # Find dataset config
    dataset_info = None
    for d in config["datasets"]:
        if d["id"] == dataset_id:
            dataset_info = d
            break
    
    if dataset_info is None:
        raise ValueError(f"Dataset {dataset_id} not found in config")
    
    # Load raw data
    logging.info(f"Loading dataset {dataset_id}")
    adata = load_matrix(dataset_info, config=config)
    
    # Compute QC metrics
    adata = compute_qc_metrics(adata)

    # Apply filters and collect mask/summary
    filtered, mask, summary = apply_filters(adata, config["qc"])

    # Persist QC artifacts
    qc_config = config["qc"].copy()
    qc_config["metrics_dir"] = config["outputs"].get("metrics_dir", "processed/metrics")
    _write_qc_tables(dataset_id, adata, mask, qc_config, summary)
    _write_qc_plots(adata[mask], dataset_id, config["qc"])

    # Save to interim
    ensure_dir(Path("data/interim"))
    outpath = f"data/interim/{dataset_id}.h5ad"
    filtered.write(outpath)

    logging.info(f"Saved QC'd dataset to {outpath}")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Run QC on datasets")
    parser.add_argument("--run", action="store_true", help="Run QC on all datasets")
    parser.add_argument("--dataset", help="Process specific dataset")
    parser.add_argument("--config", default="config/atlas.yaml", help="Config file path")
    
    args = parser.parse_args()
    
    setup_logging()
    config = load_config(args.config)
    
    if args.dataset:
        # Process single dataset
        with timer(f"QC for dataset {args.dataset}"):
            process_dataset_qc(args.dataset, config)
    elif args.run:
        # Process all datasets
        for dataset_info in config["datasets"]:
            dataset_id = dataset_info["id"]
            with timer(f"QC for dataset {dataset_id}"):
                process_dataset_qc(dataset_id, config)
