"""Doublet detection and removal for the Single-cell Immune Atlas."""

import json
import logging
from pathlib import Path
from typing import Dict, Optional, Tuple

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import scrublet as scr

from .utils import ensure_dir, load_config, setup_logging, timer


def _mask_genes(adata: ad.AnnData, enable: bool) -> ad.AnnData:
    """Optionally remove mitochondrial/ribosomal genes prior to Scrublet."""
    if not enable:
        return adata

    var_names = adata.var_names.astype(str)
    keep_mask = ~(var_names.str.startswith("MT-") | var_names.str.startswith("RPS") | var_names.str.startswith("RPL"))
    if keep_mask.all():
        return adata

    dropped_genes = var_names[~keep_mask]
    logging.info(
        "Masking %d high mitochondrial/ribosomal genes before Scrublet (examples: %s)",
        dropped_genes.size,
        ", ".join(map(str, dropped_genes[:5]))
    )
    return adata[:, keep_mask].copy()


def run_scrublet(
    adata: ad.AnnData,
    doublet_config: Dict,
    seed: int = 0,
) -> Tuple[scr.Scrublet, np.ndarray, np.ndarray]:
    """Run Scrublet doublet detection and return the fitted object and results."""
    expected_doublet_rate = doublet_config.get("expected_doublet_rate", 0.06)
    logging.info("Running Scrublet with expected doublet rate: %s", expected_doublet_rate)

    adata_masked = _mask_genes(adata, doublet_config.get("mask_high_mt_genes", False))

    scrub = scr.Scrublet(
        adata_masked.X,
        expected_doublet_rate=expected_doublet_rate,
        random_state=seed,
    )

    scrublet_kwargs = doublet_config.get("scrublet_kwargs", {})
    doublet_scores, predicted_doublets = scrub.scrub_doublets(**scrublet_kwargs)

    n_doublets = int(predicted_doublets.sum())
    doublet_rate = n_doublets / len(predicted_doublets)

    logging.info("Detected %d doublets (%.2f%% of cells)", n_doublets, 100 * doublet_rate)

    return scrub, doublet_scores, predicted_doublets


def filter_doublets(
    adata: ad.AnnData,
    scrub: scr.Scrublet,
    doublet_scores: np.ndarray,
    predicted_doublets: np.ndarray,
) -> ad.AnnData:
    """Annotate and filter predicted doublets using a fitted Scrublet instance."""
    adata = adata.copy()
    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets.astype(bool)
    threshold = float(getattr(scrub, "threshold_", np.nan))
    adata.obs["scrublet_threshold"] = threshold

    n_doublets = int(predicted_doublets.sum())
    n_cells_before = adata.n_obs

    adata_filtered = adata[~adata.obs["predicted_doublet"]].copy()
    n_cells_after = adata_filtered.n_obs

    logging.info(
        "Doublet filtering: %d → %d cells (%d doublets removed, %.1f%% retained)",
        n_cells_before,
        n_cells_after,
        n_doublets,
        100 * n_cells_after / max(n_cells_before, 1),
    )

    return adata_filtered


def _save_doublet_artifacts(
    dataset_id: str,
    adata: ad.AnnData,
    scrub: scr.Scrublet,
    doublet_scores: np.ndarray,
    predicted_doublets: np.ndarray,
    config: Dict,
) -> None:
    """Persist Scrublet diagnostics and summary files."""
    metrics_dir = Path(config["outputs"].get("metrics_dir", "processed/metrics")) / "doublets"
    ensure_dir(metrics_dir)
    figures_dir = Path(config["qc"].get("qc_plots_dir", "processed/figures/qc")).parent / "doublets"
    ensure_dir(figures_dir)

    summary = {
        "dataset_id": dataset_id,
        "n_cells": int(adata.n_obs),
        "n_doublets": int(predicted_doublets.sum()),
        "doublet_rate": float(predicted_doublets.mean()),
        "threshold": float(getattr(scrub, "threshold_", np.nan)),
        "expected_doublet_rate": float(config["doublets"].get("expected_doublet_rate", 0.06)),
    }

    summary_path = metrics_dir / f"{dataset_id}_doublet_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    logging.info("Saved doublet summary to %s", summary_path)

    fig_path = figures_dir / f"{dataset_id}_scrublet_hist.png"
    if config["doublets"].get("save_plots", True):
        fig = scrub.plot_histogram()
        fig.savefig(fig_path, dpi=300)
        plt.close(fig)
        logging.info("Saved Scrublet histogram to %s", fig_path)
    else:
        fig_path.write_text("# Plots disabled in configuration\n")


def process_dataset_doublets(dataset_id: str, config: dict) -> None:
    """Process doublet detection for a single dataset."""
    # Load QC'd data
    input_path = f"data/interim/{dataset_id}.h5ad"
    logging.info(f"Loading QC'd dataset from {input_path}")
    
    adata = ad.read_h5ad(input_path)

    seed = config.get("seed", 0)

    # Run Scrublet
    scrub, doublet_scores, predicted_doublets = run_scrublet(adata, config["doublets"], seed=seed)

    # Annotate/filter doublets
    adata_with_doublets = adata.copy()
    adata_with_doublets.obs["doublet_score"] = doublet_scores
    adata_with_doublets.obs["predicted_doublet"] = predicted_doublets.astype(bool)
    adata_with_doublets.obs["scrublet_threshold"] = float(getattr(scrub, "threshold_", np.nan))

    _save_doublet_artifacts(dataset_id, adata_with_doublets, scrub, doublet_scores, predicted_doublets, config)

    annotated_path = Path("data/interim") / f"{dataset_id}.with_doublets.h5ad"
    adata_with_doublets.write(annotated_path)
    logging.info("Saved doublet-annotated dataset to %s", annotated_path)

    adata_filtered = adata_with_doublets[~adata_with_doublets.obs["predicted_doublet"]].copy()

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
