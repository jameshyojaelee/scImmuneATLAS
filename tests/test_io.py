"""Tests for I/O functions."""

import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from src.atlas.io import download_if_needed, load_matrix


def test_download_if_needed_local_path():
    """Test that local paths are returned as-is."""
    entry = {
        "id": "test_dataset",
        "url": "/local/path/data.h5ad",
        "cancer_type": "test"
    }
    
    result = download_if_needed(entry, Path("/tmp"))
    
    assert result["url"] == "/local/path/data.h5ad"
    assert result["id"] == "test_dataset"


def test_load_matrix_metadata():
    """Test that load_matrix adds required metadata."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a dummy H5AD file
        n_obs, n_vars = 50, 30
        X = np.random.negative_binomial(5, 0.3, size=(n_obs, n_vars))
        adata = ad.AnnData(X.astype(np.float32))
        adata.var_names = [f"GENE_{i}" for i in range(n_vars)]
        adata.obs_names = [f"CELL_{i}" for i in range(n_obs)]
        
        test_path = Path(tmpdir) / "test.h5ad"
        adata.write(test_path)
        
        # Test loading with metadata
        entry = {
            "id": "test_dataset",
            "url": str(test_path),
            "cancer_type": "melanoma",
            "platform": "10x"
        }
        
        loaded_adata = load_matrix(entry)
        
        # Check metadata was added
        assert "dataset_id" in loaded_adata.obs.columns
        assert "cancer_type" in loaded_adata.obs.columns
        assert "platform" in loaded_adata.obs.columns
        
        assert loaded_adata.obs["dataset_id"].iloc[0] == "test_dataset"
        assert loaded_adata.obs["cancer_type"].iloc[0] == "melanoma"
        assert loaded_adata.obs["platform"].iloc[0] == "10x"


def test_load_matrix_unique_names():
    """Test that var_names and obs_names are made unique."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create data with duplicate names
        n_obs, n_vars = 20, 10
        X = np.random.randn(n_obs, n_vars)
        adata = ad.AnnData(X)
        
        # Add duplicate gene names
        adata.var_names = ["GENE_A"] * 5 + ["GENE_B"] * 5
        adata.obs_names = ["CELL_1"] * 10 + ["CELL_2"] * 10
        
        test_path = Path(tmpdir) / "test.h5ad"
        adata.write(test_path)
        
        entry = {
            "id": "test_dataset",
            "url": str(test_path),
            "cancer_type": "test"
        }
        
        loaded_adata = load_matrix(entry)
        
        # Names should be unique after loading
        assert len(set(loaded_adata.var_names)) == len(loaded_adata.var_names)
        assert len(set(loaded_adata.obs_names)) == len(loaded_adata.obs_names)
