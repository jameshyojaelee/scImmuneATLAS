"""Doublet detection and removal for the Single-cell Immune Atlas."""

import logging
from pathlib import Path
from typing import Optional

import anndata as ad
import numpy as np
import scrublet as scr

from .utils import load_config, setup_logging, timer


def run_scrublet(adata: ad.AnnData, expected_doublet_rate: float = 0.06) -> np.ndarray:
    """Run Scrublet doublet detection."""
    logging.info(f"Running Scrublet with expected doublet rate: {expected_doublet_rate}")
    
    # Initialize Scrublet
    scrub = scr.Scrublet(
        adata.X, 
        expected_doublet_rate=expected_doublet_rate,
        random_state=0
    )
    
    # Run doublet detection
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=2, 
        min_cells=3, 
        min_gene_variability_pctl=85, 
        n_prin_comps=30
    )
    
    n_doublets = predicted_doublets.sum()
    doublet_rate = n_doublets / len(predicted_doublets)
    
    logging.info(f"Detected {n_doublets} doublets ({doublet_rate:.2%} of cells)")
    
    return doublet_scores


def filter_doublets(
    adata: ad.AnnData, 
    doublet_scores: np.ndarray, 
    threshold: str = "auto"
) -> ad.AnnData:
    """Filter out predicted doublets."""
    # Add doublet scores to obs
    adata.obs["doublet_score"] = doublet_scores
    
    if threshold == "auto":
        # Use Scrublet's automatic threshold
        scrub = scr.Scrublet(adata.X, random_state=0)
        _, predicted_doublets = scrub.scrub_doublets()
        adata.obs["predicted_doublet"] = predicted_doublets
    else:
        # Use manual threshold
        adata.obs["predicted_doublet"] = doublet_scores > threshold
    
    n_doublets = adata.obs["predicted_doublet"].sum()
    n_cells_before = adata.n_obs
    
    # Filter out doublets
    adata_filtered = adata[~adata.obs["predicted_doublet"]].copy()
    
    n_cells_after = adata_filtered.n_obs
    
    logging.info(f"Doublet filtering: {n_cells_before} → {n_cells_after} cells "
                f"({n_doublets} doublets removed, {100 * n_cells_after / n_cells_before:.1f}% retained)")
    
    return adata_filtered


def process_dataset_doublets(dataset_id: str, config: dict) -> None:
    """Process doublet detection for a single dataset."""
    # Load QC'd data
    input_path = f"data/interim/{dataset_id}.h5ad"
    logging.info(f"Loading QC'd dataset from {input_path}")
    
    adata = ad.read_h5ad(input_path)
    
    # Run Scrublet
    doublet_scores = run_scrublet(adata, config["doublets"]["expected_doublet_rate"])
    
    # Filter doublets
    adata_filtered = filter_doublets(adata, doublet_scores)
    
    # Save filtered data
    output_path = f"data/interim/{dataset_id}.doublet_filtered.h5ad"
    adata_filtered.write(output_path)
    
    logging.info(f"Saved doublet-filtered dataset to {output_path}")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Run doublet detection")
    parser.add_argument("--run", action="store_true", help="Run doublet detection on all datasets")
    parser.add_argument("--dataset", help="Process specific dataset")
    parser.add_argument("--config", default="config/atlas.yaml", help="Config file path")
    
    args = parser.parse_args()
    
    setup_logging()
    config = load_config(args.config)
    
    if args.dataset:
        # Process single dataset
        with timer(f"Doublet detection for dataset {args.dataset}"):
            process_dataset_doublets(args.dataset, config)
    elif args.run:
        # Process all datasets
        for dataset_info in config["datasets"]:
            dataset_id = dataset_info["id"]
            with timer(f"Doublet detection for dataset {dataset_id}"):
                process_dataset_doublets(dataset_id, config)
