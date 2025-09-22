"""Smoke tests for the Single-cell Immune Atlas."""

import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pytest

from src.atlas import annotate, doublets, export, integration, io, qc, utils, viz


def test_imports():
    """Test that all modules can be imported."""
    # All modules should import without error
    assert hasattr(utils, "load_config")
    assert hasattr(io, "load_matrix")
    assert hasattr(qc, "compute_qc_metrics")
    assert hasattr(doublets, "run_scrublet")
    assert hasattr(integration, "integrate_scvi")
    assert hasattr(annotate, "load_immune_markers")
    assert hasattr(viz, "umap_by")
    assert hasattr(export, "write_cellxgene")


def test_config_loading():
    """Test configuration loading."""
    config = utils.load_config("config/atlas.yaml")
    
    assert "project_name" in config
    assert "datasets" in config
    assert "qc" in config
    assert "integration" in config
    
    # Check QC parameters
    assert "min_genes" in config["qc"]
    assert "max_genes" in config["qc"]
    assert "max_mt_pct" in config["qc"]


def test_markers_loading():
    """Test immune markers loading."""
    markers_df = annotate.load_immune_markers()
    
    assert len(markers_df) > 0
    assert "cell_type" in markers_df.columns
    assert "gene" in markers_df.columns
    assert "direction" in markers_df.columns
    
    # Check for key cell types
    cell_types = set(markers_df["cell_type"])
    expected_types = {"CD8_T", "CD4_T", "NK", "B_cell", "Mono"}
    assert len(expected_types & cell_types) > 0


def test_dummy_adata_creation():
    """Test creating a dummy AnnData object."""
    # Create synthetic data
    n_obs, n_vars = 100, 50
    X = np.random.negative_binomial(5, 0.3, size=(n_obs, n_vars))
    adata = ad.AnnData(X.astype(np.float32))
    
    # Add gene names with some MT genes
    gene_names = [f"GENE_{i}" for i in range(45)]
    gene_names.extend([f"MT-{i}" for i in range(5)])
    adata.var_names = gene_names
    
    # Add cell metadata
    adata.obs["dataset_id"] = "test_dataset"
    adata.obs["cancer_type"] = "test_cancer"
    
    assert adata.n_obs == n_obs
    assert adata.n_vars == n_vars
    assert "dataset_id" in adata.obs.columns


def test_qc_metrics():
    """Test QC metrics computation."""
    # Create dummy data
    n_obs, n_vars = 50, 30
    X = np.random.negative_binomial(5, 0.3, size=(n_obs, n_vars))
    adata = ad.AnnData(X.astype(np.float32))
    
    # Add gene names with MT genes
    gene_names = [f"GENE_{i}" for i in range(25)]
    gene_names.extend([f"MT-{i}" for i in range(5)])
    adata.var_names = gene_names
    
    # Compute QC metrics
    adata_qc = qc.compute_qc_metrics(adata)
    
    assert "pct_counts_mt" in adata_qc.obs.columns
    assert "n_genes_by_counts" in adata_qc.obs.columns
    assert "total_counts" in adata_qc.obs.columns
    
    # Check that MT percentage is reasonable
    assert adata_qc.obs["pct_counts_mt"].max() <= 100
    assert adata_qc.obs["pct_counts_mt"].min() >= 0


def test_simple_viz():
    """Test simple visualization functions."""
    # Create dummy data with UMAP coordinates
    n_obs = 100
    adata = ad.AnnData(np.random.randn(n_obs, 10))
    adata.obs["cell_type"] = np.random.choice(["T_cell", "B_cell", "NK"], n_obs)
    adata.obsm["X_umap"] = np.random.randn(n_obs, 2)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        outpath = Path(tmpdir) / "test_umap.png"
        
        # This should not raise an error
        viz.umap_by(adata, "cell_type", str(outpath))
        
        # Check that file was created
        assert outpath.exists()


def test_seed_setting():
    """Test seed setting for reproducibility."""
    utils.set_seed(42)
    
    # Generate some random numbers
    rand1 = np.random.randn(10)
    
    # Reset seed
    utils.set_seed(42)
    rand2 = np.random.randn(10)
    
    # Should be identical
    np.testing.assert_array_equal(rand1, rand2)
