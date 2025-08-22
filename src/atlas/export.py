"""Export functions for the Single-cell Immune Atlas."""

import logging
from pathlib import Path
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd

from .utils import ensure_dir, load_config, setup_logging, timer


def sanitize_for_cellxgene(adata: ad.AnnData) -> ad.AnnData:
    """Sanitize AnnData object for cellxgene compatibility."""
    # Make a copy
    adata = adata.copy()
    
    # Ensure all obs columns are strings or categoricals
    for col in adata.obs.columns:
        if adata.obs[col].dtype == "object":
            # Convert to categorical for efficiency
            adata.obs[col] = adata.obs[col].astype("category")
        elif pd.api.types.is_numeric_dtype(adata.obs[col]):
            # Keep numeric columns as-is, but ensure no NaNs
            adata.obs[col] = adata.obs[col].fillna(0)
    
    # Ensure var_names are strings
    adata.var_names = adata.var_names.astype(str)
    adata.obs_names = adata.obs_names.astype(str)
    
    # Remove any sparse categories with NaN
    for col in adata.obs.columns:
        if pd.api.types.is_categorical_dtype(adata.obs[col]):
            # Remove unused categories
            adata.obs[col] = adata.obs[col].cat.remove_unused_categories()
    
    # Ensure X is not too sparse for cellxgene
    if hasattr(adata.X, "toarray"):
        # Keep as sparse if reasonable sparsity
        sparsity = 1.0 - (adata.X.nnz / (adata.n_obs * adata.n_vars))
        logging.info(f"Data sparsity: {sparsity:.2%}")
    
    # Add required metadata for cellxgene
    if "cell_type" not in adata.obs.columns:
        adata.obs["cell_type"] = "Unknown"
    
    logging.info("Sanitized data for cellxgene compatibility")
    
    return adata


def write_cellxgene(adata: ad.AnnData, outdir: str) -> None:
    """Write cellxgene-compatible H5AD file."""
    outdir = Path(outdir)
    ensure_dir(outdir)
    
    # Sanitize data
    adata_clean = sanitize_for_cellxgene(adata)
    
    # Write with compression
    outpath = outdir / "atlas.h5ad"
    adata_clean.write(outpath, compression="gzip")
    
    logging.info(f"Exported cellxgene-compatible atlas to {outpath}")
    
    # Write metadata summary
    metadata_path = outdir / "metadata_summary.txt"
    with open(metadata_path, "w") as f:
        f.write("Single-cell Immune Atlas - Cellxgene Export\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Cells: {adata_clean.n_obs:,}\n")
        f.write(f"Genes: {adata_clean.n_vars:,}\n\n")
        
        f.write("Cell type distribution:\n")
        type_counts = adata_clean.obs["cell_type"].value_counts()
        for cell_type, count in type_counts.items():
            pct = 100 * count / adata_clean.n_obs
            f.write(f"  {cell_type}: {count:,} ({pct:.1f}%)\n")
        
        if "dataset_id" in adata_clean.obs.columns:
            f.write("\nDataset distribution:\n")
            dataset_counts = adata_clean.obs["dataset_id"].value_counts()
            for dataset, count in dataset_counts.items():
                pct = 100 * count / adata_clean.n_obs
                f.write(f"  {dataset}: {count:,} ({pct:.1f}%)\n")
        
        if "cancer_type" in adata_clean.obs.columns:
            f.write("\nCancer type distribution:\n")
            cancer_counts = adata_clean.obs["cancer_type"].value_counts()
            for cancer, count in cancer_counts.items():
                pct = 100 * count / adata_clean.n_obs
                f.write(f"  {cancer}: {count:,} ({pct:.1f}%)\n")
    
    logging.info(f"Wrote metadata summary to {metadata_path}")


def write_seurat(adata: ad.AnnData, outpath: str) -> None:
    """Write Seurat-compatible object (requires R and SeuratDisk)."""
    try:
        import anndata2ri
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        
        # Activate conversions
        anndata2ri.activate()
        pandas2ri.activate()
        
        # Convert to Seurat object in R
        ro.r(f'''
        library(SeuratDisk)
        library(Seurat)
        
        # This would require the AnnData object to be passed to R
        # For now, we'll save as h5ad and convert in R separately
        ''')
        
        logging.info("Seurat export requires manual conversion in R")
        logging.info(f"Use: Convert('{outpath}.h5ad', dest='h5seurat', overwrite=TRUE)")
        
    except ImportError:
        logging.warning("R/rpy2/anndata2ri not available for Seurat export")
        logging.info("To enable Seurat export, install: pip install rpy2 anndata2ri")


def run_export(config: dict) -> None:
    """Run export based on config settings."""
    # Load annotated data
    input_path = "processed/integrated_annotated.h5ad"
    logging.info(f"Loading annotated data from {input_path}")
    
    adata = ad.read_h5ad(input_path)
    
    # Export for cellxgene
    cellxgene_dir = config["outputs"]["cellxgene_dir"]
    write_cellxgene(adata, cellxgene_dir)
    
    # Optionally export for Seurat
    # write_seurat(adata, "processed/atlas")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Export atlas for viewers")
    parser.add_argument("--cellxgene", action="store_true", help="Export for cellxgene")
    parser.add_argument("--seurat", action="store_true", help="Export for Seurat")
    parser.add_argument("--config", default="config/atlas.yaml", help="Config file path")
    
    args = parser.parse_args()
    
    setup_logging()
    config = load_config(args.config)
    
    if args.cellxgene:
        with timer("Cellxgene export"):
            run_export(config)
    
    if args.seurat:
        with timer("Seurat export"):
            input_path = "processed/integrated_annotated.h5ad"
            adata = ad.read_h5ad(input_path)
            write_seurat(adata, "processed/atlas")
