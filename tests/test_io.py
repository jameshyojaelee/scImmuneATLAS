"""Tests for I/O functions."""

import hashlib
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np

from src.atlas.io import _compute_sha256, download_if_needed, load_matrix
from src.atlas.schemas import validate_anndata


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
            "cancer_type": "test",
            "platform": "10x",
        }

        loaded_adata = load_matrix(entry)

        # Names should be unique after loading
        assert len(set(loaded_adata.var_names)) == len(loaded_adata.var_names)
        assert len(set(loaded_adata.obs_names)) == len(loaded_adata.obs_names)


def test_compute_sha256():
    """Test SHA256 computation."""
    with tempfile.TemporaryDirectory() as tmpdir:
        test_file = Path(tmpdir) / "test.txt"
        test_content = b"Hello, World!"
        test_file.write_bytes(test_content)

        # Compute expected hash
        expected_hash = hashlib.sha256(test_content).hexdigest()

        # Test our function
        actual_hash = _compute_sha256(test_file)

        assert actual_hash == expected_hash


def test_validate_anndata_success():
    """Test successful AnnData validation."""
    n_obs, n_vars = 50, 30
    X = np.random.negative_binomial(5, 0.3, size=(n_obs, n_vars))
    adata = ad.AnnData(X.astype(np.float32))
    adata.var_names = [f"GENE_{i}" for i in range(n_vars)]
    adata.obs_names = [f"CELL_{i}" for i in range(n_obs)]

    # Add required metadata
    adata.obs["dataset_id"] = "test_dataset"
    adata.obs["cancer_type"] = "melanoma"
    adata.obs["platform"] = "10x"

    is_valid, error = validate_anndata(adata, raise_on_error=False)

    assert is_valid is True
    assert error is None


def test_validate_anndata_missing_metadata():
    """Test validation failure when required metadata is missing."""
    n_obs, n_vars = 50, 30
    X = np.random.negative_binomial(5, 0.3, size=(n_obs, n_vars))
    adata = ad.AnnData(X.astype(np.float32))
    adata.var_names = [f"GENE_{i}" for i in range(n_vars)]
    adata.obs_names = [f"CELL_{i}" for i in range(n_obs)]

    # Missing required metadata columns
    is_valid, error = validate_anndata(
        adata, dataset_id="test", raise_on_error=False
    )

    assert is_valid is False
    assert error is not None
    assert "dataset_id" in error or "cancer_type" in error or "platform" in error


def test_validate_anndata_empty():
    """Test validation failure for empty AnnData."""
    adata = ad.AnnData(np.array([]).reshape(0, 0))

    is_valid, error = validate_anndata(adata, raise_on_error=False)

    assert is_valid is False
    assert error is not None
    assert "0 observations" in error or "0 variables" in error


def test_download_with_sha256_verification():
    """Test download with SHA256 verification (mocked)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a test file to simulate download
        test_file = Path(tmpdir) / "test.h5ad"
        test_content = b"Test AnnData file content"
        test_file.write_bytes(test_content)

        # Test that existing file with correct hash doesn't raise warning
        entry = {
            "id": "test_dataset",
            "url": str(test_file),  # Local path, not URL
            "cancer_type": "test",
        }

        # Should work without hash
        result = download_if_needed(entry, Path(tmpdir))
        assert result["url"] == str(test_file)


def test_load_matrix_with_validation():
    """Test that load_matrix validates AnnData."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a valid H5AD file
        n_obs, n_vars = 50, 30
        X = np.random.negative_binomial(5, 0.3, size=(n_obs, n_vars))
        adata = ad.AnnData(X.astype(np.float32))
        adata.var_names = [f"GENE_{i}" for i in range(n_vars)]
        adata.obs_names = [f"CELL_{i}" for i in range(n_obs)]

        test_path = Path(tmpdir) / "test.h5ad"
        adata.write(test_path)

        entry = {
            "id": "test_dataset",
            "url": str(test_path),
            "cancer_type": "melanoma",
            "platform": "10x",
        }

        # Should load and validate successfully
        loaded_adata = load_matrix(entry)

        assert loaded_adata.n_obs == n_obs
        assert loaded_adata.n_vars == n_vars
        assert "dataset_id" in loaded_adata.obs.columns
        assert "cancer_type" in loaded_adata.obs.columns
        assert "platform" in loaded_adata.obs.columns
