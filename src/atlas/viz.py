"""Visualization functions for the Single-cell Immune Atlas."""

import logging
from pathlib import Path
from typing import List, Optional

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
import seaborn as sns

from .utils import ensure_dir, load_config, setup_logging, timer


def umap_by(adata: ad.AnnData, key: str, outpath: str, figsize: tuple = (8, 6)) -> None:
    """Create UMAP plot colored by a metadata key."""
    ensure_dir(Path(outpath).parent)
    
    # Set up the plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot UMAP
    sc.pl.umap(adata, color=key, ax=ax, show=False, frameon=False)
    
    # Clean up the plot
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(f"UMAP colored by {key}")
    
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()
    
    logging.info(f"Saved UMAP plot to {outpath}")


def stacked_bar(
    adata: ad.AnnData, 
    groupby: List[str], 
    normalize: bool = True, 
    outpath: str,
    figsize: tuple = (10, 6)
) -> None:
    """Create stacked bar plot showing proportions."""
    ensure_dir(Path(outpath).parent)
    
    # Create crosstab
    if len(groupby) == 2:
        crosstab = pd.crosstab(adata.obs[groupby[0]], adata.obs[groupby[1]])
    else:
        raise ValueError("groupby must contain exactly 2 elements")
    
    # Normalize if requested
    if normalize:
        crosstab = crosstab.div(crosstab.sum(axis=1), axis=0)
    
    # Create plot
    fig, ax = plt.subplots(figsize=figsize)
    crosstab.plot(kind="bar", stacked=True, ax=ax, colormap="tab20")
    
    ax.set_xlabel(groupby[0])
    ax.set_ylabel("Proportion" if normalize else "Count")
    ax.set_title(f"{groupby[1]} by {groupby[0]}")
    ax.legend(title=groupby[1], bbox_to_anchor=(1.05, 1), loc="upper left")
    
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()
    
    logging.info(f"Saved stacked bar plot to {outpath}")


def marker_heatmap(
    adata: ad.AnnData, 
    cell_types: List[str], 
    genes: List[str], 
    outpath: str,
    figsize: tuple = (12, 8)
) -> None:
    """Create marker heatmap showing mean expression per cell type."""
    ensure_dir(Path(outpath).parent)
    
    # Filter to available genes
    available_genes = [g for g in genes if g in adata.var_names]
    if len(available_genes) == 0:
        logging.warning("No marker genes found in dataset")
        return
    
    # Filter to available cell types
    available_types = [ct for ct in cell_types if ct in adata.obs["cell_type"].values]
    if len(available_types) == 0:
        logging.warning("No specified cell types found in dataset")
        return
    
    # Compute mean expression per cell type
    expr_data = []
    for cell_type in available_types:
        mask = adata.obs["cell_type"] == cell_type
        if mask.sum() == 0:
            continue
        
        subset = adata[mask, available_genes]
        if hasattr(subset.X, "toarray"):
            mean_expr = np.mean(subset.X.toarray(), axis=0)
        else:
            mean_expr = np.mean(subset.X, axis=0)
        
        expr_data.append(mean_expr)
    
    if len(expr_data) == 0:
        logging.warning("No data to plot")
        return
    
    # Create DataFrame
    expr_df = pd.DataFrame(
        expr_data, 
        index=available_types[:len(expr_data)], 
        columns=available_genes
    )
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(expr_df, annot=True, fmt=".2f", cmap="viridis", ax=ax)
    
    ax.set_title("Marker Gene Expression by Cell Type")
    ax.set_xlabel("Genes")
    ax.set_ylabel("Cell Types")
    
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()
    
    logging.info(f"Saved marker heatmap to {outpath}")


def generate_all_figures(config: dict) -> None:
    """Generate all standard figures."""
    # Load annotated data
    input_path = "processed/integrated_annotated.h5ad"
    logging.info(f"Loading annotated data from {input_path}")
    
    adata = ad.read_h5ad(input_path)
    
    # Create output directory
    figures_dir = Path(config["outputs"]["figures_dir"])
    ensure_dir(figures_dir)
    
    # Generate UMAP plots
    umap_by(adata, "cell_type", figures_dir / "umap_by_cell_type.png")
    umap_by(adata, "dataset_id", figures_dir / "umap_by_dataset.png")
    
    if "cancer_type" in adata.obs.columns:
        umap_by(adata, "cancer_type", figures_dir / "umap_by_cancer_type.png")
    
    # Generate proportion plots
    if "cancer_type" in adata.obs.columns:
        stacked_bar(
            adata, 
            ["cancer_type", "cell_type"], 
            True, 
            figures_dir / "proportions_by_cancer_type.png"
        )
    
    stacked_bar(
        adata, 
        ["dataset_id", "cell_type"], 
        True, 
        figures_dir / "proportions_by_dataset.png"
    )
    
    # Generate marker heatmap
    key_markers = [
        "CD8A", "CD4", "FOXP3", "NKG7", "MS4A1", "MZB1", 
        "LST1", "CLEC9A", "FCER1A"
    ]
    key_types = ["CD8_T", "CD4_T", "Treg", "NK", "B_cell", "Plasma", "Mono", "DC_cDC1", "DC_cDC2"]
    
    marker_heatmap(
        adata, 
        key_types, 
        key_markers, 
        figures_dir / "marker_heatmap.png"
    )


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate visualizations")
    parser.add_argument("--all", action="store_true", help="Generate all figures")
    parser.add_argument("--config", default="config/atlas.yaml", help="Config file path")
    
    args = parser.parse_args()
    
    setup_logging()
    config = load_config(args.config)
    
    if args.all:
        with timer("Generating all figures"):
            generate_all_figures(config)
