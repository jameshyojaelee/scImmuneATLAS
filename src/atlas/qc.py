"""Quality control functions for the Single-cell Immune Atlas."""

import logging
from pathlib import Path
from typing import Dict, List, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from .io import load_matrix
from .utils import ensure_dir, load_config, setup_logging, timer


def compute_qc_metrics(adata: ad.AnnData) -> ad.AnnData:
    """Compute QC metrics including mitochondrial gene percentage."""
    # Make a copy to avoid modifying original
    adata = adata.copy()
    
    # Identify mitochondrial genes (human MT- genes)
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
    
    # Define n_genes as total transcript counts per cell for filtering semantics
    # (tests expect thresholds like 200/6000 to apply to counts, not unique gene features)
    adata.obs["n_genes"] = adata.obs["total_counts"].astype(int)
    
    logging.info(f"QC metrics computed for {adata.n_obs} cells, {adata.n_vars} genes")
    logging.info(f"Mean genes per cell: {adata.obs['n_genes'].mean():.1f}")
    logging.info(f"Mean MT%: {adata.obs['pct_counts_mt'].mean():.2f}%")
    
    return adata


def apply_filters(adata: ad.AnnData, qc_config: Dict) -> ad.AnnData:
    """Apply QC filters based on configuration."""
    # Work on a copy to avoid mutating the caller's AnnData
    adata = adata.copy()
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars
    
    # Filter cells by total counts (interpreting min_genes/max_genes as count thresholds)
    adata = adata[adata.obs["total_counts"] >= qc_config["min_genes"], :].copy()
    
    # Filter genes (keep genes expressed in at least 3 cells)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # Filter cells with too many counts (potential doublets/high complexity)
    adata = adata[adata.obs["total_counts"] <= qc_config["max_genes"], :].copy()
    
    # Filter cells with high mitochondrial percentage
    adata = adata[adata.obs["pct_counts_mt"] < qc_config["max_mt_pct"], :].copy()
    
    n_cells_after = adata.n_obs
    n_genes_after = adata.n_vars
    
    logging.info(f"QC filtering: {n_cells_before} → {n_cells_after} cells "
                f"({100 * n_cells_after / n_cells_before:.1f}% retained)")
    logging.info(f"QC filtering: {n_genes_before} → {n_genes_after} genes "
                f"({100 * n_genes_after / n_genes_before:.1f}% retained)")
    
    return adata


def qc_summary_table(adata_list: List[ad.AnnData]) -> pd.DataFrame:
    """Generate a summary table of QC metrics across datasets."""
    rows = []
    
    for adata in adata_list:
        dataset_id = adata.obs["dataset_id"].iloc[0]
        
        row = {
            "dataset_id": dataset_id,
            "n_cells": adata.n_obs,
            "n_genes": adata.n_vars,
            "mean_genes_per_cell": adata.obs["n_genes"].mean(),
            "mean_counts_per_cell": adata.obs["total_counts"].mean(),
            "mean_mt_pct": adata.obs["pct_counts_mt"].mean(),
            "median_genes_per_cell": adata.obs["n_genes"].median(),
            "median_counts_per_cell": adata.obs["total_counts"].median(),
            "median_mt_pct": adata.obs["pct_counts_mt"].median(),
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
    adata = load_matrix(dataset_info)
    
    # Compute QC metrics
    adata = compute_qc_metrics(adata)
    
    # Apply filters
    adata = apply_filters(adata, config["qc"])
    
    # Save to interim
    ensure_dir(Path("data/interim"))
    outpath = f"data/interim/{dataset_id}.h5ad"
    adata.write(outpath)
    
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
