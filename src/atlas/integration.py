"""Integration functions for the Single-cell Immune Atlas."""

import logging
from pathlib import Path
from typing import List, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.preprocessing import LabelEncoder

from .utils import load_config, set_seed, setup_logging, timer


def integrate_scvi(
    adatas: List[ad.AnnData], 
    batch_key: str, 
    latent_dim: int = 30,
    max_epochs: int = 200
) -> ad.AnnData:
    """Integrate datasets using scVI."""
    try:
        import scvi
    except ImportError:
        raise ImportError("scvi-tools is required for scVI integration")
    
    logging.info(f"Starting scVI integration with {len(adatas)} datasets")
    
    # Concatenate datasets
    adata = ad.concat(adatas, join="outer", index_unique="-")
    
    # Make var_names unique and fill missing values
    adata.var_names_unique()
    
    # Fill missing values with 0
    if hasattr(adata.X, "toarray"):
        X_dense = adata.X.toarray()
    else:
        X_dense = adata.X
    
    X_dense = np.nan_to_num(X_dense, nan=0.0)
    adata.X = X_dense
    
    # Setup scVI
    logging.info("Setting up scVI model")
    scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key)
    
    # Create and train model
    model = scvi.model.SCVI(adata, n_latent=latent_dim)
    model.train(max_epochs=max_epochs, early_stopping=True, use_gpu=False)
    
    # Get latent representation
    logging.info("Computing latent representation")
    adata.obsm["X_scvi"] = model.get_latent_representation()
    
    # Compute neighbors and UMAP on scVI representation
    sc.pp.neighbors(adata, use_rep="X_scvi", n_neighbors=15, random_state=0)
    sc.tl.umap(adata, random_state=0, min_dist=0.4, spread=1.0)
    
    logging.info(f"scVI integration completed: {adata.n_obs} cells, {adata.n_vars} genes")
    
    return adata


def integrate_harmony(adata: ad.AnnData, batch_key: str) -> ad.AnnData:
    """Integrate datasets using Harmony."""
    try:
        import harmonypy as hm
    except ImportError:
        raise ImportError("harmonypy is required for Harmony integration")
    
    logging.info("Starting Harmony integration")
    
    # Normalize and log-transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Find highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    
    # Scale data
    sc.pp.scale(adata, max_value=10)
    
    # PCA
    sc.tl.pca(adata, svd_solver="arpack", random_state=0)
    
    # Run Harmony
    logging.info("Running Harmony batch correction")
    harmony_out = hm.run_harmony(
        adata.obsm["X_pca"], 
        adata.obs, 
        batch_key,
        random_state=0
    )
    adata.obsm["X_harmony"] = harmony_out.Z_corr.T
    
    # Compute neighbors and UMAP on Harmony representation
    sc.pp.neighbors(adata, use_rep="X_harmony", n_neighbors=15, random_state=0)
    sc.tl.umap(adata, random_state=0, min_dist=0.4, spread=1.0)
    
    logging.info(f"Harmony integration completed: {adata.n_obs} cells, {adata.n_vars} genes")
    
    return adata


def run_integration(config: dict) -> None:
    """Run integration based on config settings."""
    # Load doublet-filtered datasets
    adatas = []
    for dataset_info in config["datasets"]:
        dataset_id = dataset_info["id"]
        input_path = f"data/interim/{dataset_id}.doublet_filtered.h5ad"
        
        logging.info(f"Loading {input_path}")
        adata = ad.read_h5ad(input_path)
        adatas.append(adata)
    
    # Set seed for reproducibility
    set_seed(config["seed"])
    
    # Choose integration method
    method = config["integration"]["method"]
    batch_key = config["integration"]["batch_key"]
    latent_dim = config["integration"]["latent_dim"]
    
    if method == "scvi":
        adata_integrated = integrate_scvi(adatas, batch_key, latent_dim)
    elif method == "harmony":
        # For Harmony, concatenate first
        adata_concat = ad.concat(adatas, join="outer", index_unique="-")
        adata_integrated = integrate_harmony(adata_concat, batch_key)
    else:
        raise ValueError(f"Unknown integration method: {method}")
    
    # Save integrated data
    output_path = "processed/integrated_atlas.h5ad"
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    adata_integrated.write(output_path)
    
    logging.info(f"Saved integrated atlas to {output_path}")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Run dataset integration")
    parser.add_argument("--method", choices=["scvi", "harmony"], help="Integration method")
    parser.add_argument("--config", default="config/atlas.yaml", help="Config file path")
    
    args = parser.parse_args()
    
    setup_logging()
    config = load_config(args.config)
    
    # Override method if specified
    if args.method:
        config["integration"]["method"] = args.method
    
    with timer(f"Integration using {config['integration']['method']}"):
        run_integration(config)
