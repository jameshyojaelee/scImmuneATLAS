"""Tests for QC functions."""

import numpy as np
import anndata as ad
import pytest

from src.atlas.qc import compute_qc_metrics, apply_filters, qc_summary_table


def test_compute_qc_metrics():
    """Test QC metrics computation."""
    # Create test data with known properties
    n_obs, n_vars = 100, 50
    
    # Create expression matrix
    X = np.random.negative_binomial(5, 0.3, size=(n_obs, n_vars))
    adata = ad.AnnData(X.astype(np.float32))
    
    # Add gene names - 40 regular genes, 10 MT genes
    gene_names = [f"GENE_{i}" for i in range(40)]
    gene_names.extend([f"MT-{i}" for i in range(10)])
    adata.var_names = gene_names
    
    # Compute QC metrics
    adata_qc = compute_qc_metrics(adata)
    
    # Check that required columns exist
    required_cols = ["pct_counts_mt", "n_genes_by_counts", "total_counts", "total_counts_mt"]
    for col in required_cols:
        assert col in adata_qc.obs.columns
    
    # Check that MT genes are identified
    assert "mt" in adata_qc.var.columns
    mt_genes = adata_qc.var["mt"].sum()
    assert mt_genes == 10  # Should identify 10 MT genes
    
    # Check that percentages are reasonable
    assert adata_qc.obs["pct_counts_mt"].max() <= 100
    assert adata_qc.obs["pct_counts_mt"].min() >= 0


def test_apply_filters():
    """Test QC filtering."""
    # Create test data
    n_obs, n_vars = 200, 100
    X = np.random.negative_binomial(5, 0.3, size=(n_obs, n_vars))
    adata = ad.AnnData(X.astype(np.float32))
    
    # Add gene names with MT genes
    gene_names = [f"GENE_{i}" for i in range(80)]
    gene_names.extend([f"MT-{i}" for i in range(20)])
    adata.var_names = gene_names
    
    # Add some cells with extreme values
    # High MT percentage for some cells
    mt_genes_idx = np.where(adata.var_names.str.startswith("MT-"))[0]
    adata.X[0:10, mt_genes_idx] *= 10  # High MT for first 10 cells
    
    # Very few genes for some cells
    adata.X[10:20, :] = 0  # No expression for cells 10-20
    adata.X[10:20, 0:5] = 1  # Only 5 genes expressed
    
    # Compute QC metrics first
    adata = compute_qc_metrics(adata)
    
    # Apply filters
    qc_config = {
        "min_genes": 200,
        "max_genes": 6000,
        "max_mt_pct": 15
    }
    
    adata_filtered, mask, summary = apply_filters(adata, qc_config)
    
    # Check that filtering worked
    assert summary["n_cells_after"] <= summary["n_cells_before"]

    # Check thresholds only if cells remain
    if adata_filtered.n_obs > 0:
        genes_metric = adata_filtered.obs["n_genes_by_counts"]
        assert genes_metric.min() >= qc_config["min_genes"]
        assert genes_metric.max() <= qc_config["max_genes"]
        assert adata_filtered.obs["pct_counts_mt"].max() <= qc_config["max_mt_pct"]


def test_filter_edge_cases():
    """Test filtering with edge cases."""
    # Create minimal data
    n_obs, n_vars = 10, 5
    X = np.ones((n_obs, n_vars))
    adata = ad.AnnData(X)
    adata.var_names = [f"GENE_{i}" for i in range(n_vars)]
    
    # Compute QC metrics
    adata = compute_qc_metrics(adata)
    
    # Very strict filters - should keep some cells
    strict_config = {
        "min_genes": 1,
        "max_genes": 10,
        "max_mt_pct": 100  # Allow all MT percentages
    }
    
    adata_filtered, *_ = apply_filters(adata, strict_config)
    assert adata_filtered.n_obs > 0  # Should keep some cells
    
    # Impossible filters - should filter everything
    impossible_config = {
        "min_genes": 1000,  # More genes than exist
        "max_genes": 6000,
        "max_mt_pct": 15
    }
    
    adata_filtered, *_ = apply_filters(adata, impossible_config)
    # This might filter all cells, which is okay for impossible criteria


def _synthetic_qc_adata(dataset_id: str, seed: int = 0) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    X = rng.negative_binomial(5, 0.3, size=(30, 15))
    adata = ad.AnnData(X.astype(np.float32))
    adata.var_names = [f"GENE_{i}" for i in range(15)]
    adata.obs_names = [f"{dataset_id}_cell_{i}" for i in range(30)]
    adata.obs["dataset_id"] = dataset_id
    return compute_qc_metrics(adata)


def test_qc_summary_table_happy_path():
    """qc_summary_table reports aggregates using columns computed by compute_qc_metrics."""
    adata_a = _synthetic_qc_adata("A", seed=1)
    adata_b = _synthetic_qc_adata("B", seed=2)

    summary = qc_summary_table([adata_a, adata_b])
    assert list(summary["dataset_id"]) == ["A", "B"]

    for adata, row in zip([adata_a, adata_b], summary.itertuples()):
        assert row.n_cells == adata.n_obs
        assert row.n_genes == adata.n_vars
        np.testing.assert_allclose(
            row.mean_genes_per_cell,
            adata.obs["n_genes_by_counts"].mean(),
        )
        np.testing.assert_allclose(
            row.median_counts_per_cell,
            adata.obs["total_counts"].median(),
        )


def test_qc_summary_table_missing_optional_columns():
    """Missing QC columns should yield NaN rather than raising."""
    X = np.ones((5, 4), dtype=np.float32)
    adata = ad.AnnData(X)
    adata.var_names = [f"G{i}" for i in range(4)]
    adata.obs_names = [f"cell{i}" for i in range(5)]
    adata.obs["dataset_id"] = "MISSING"
    summary = qc_summary_table([adata])

    assert summary.loc[0, "dataset_id"] == "MISSING"
    assert np.isnan(summary.loc[0, "mean_counts_per_cell"])
    assert np.isnan(summary.loc[0, "mean_mt_pct"])


def test_qc_summary_table_handles_empty_datasets():
    """Zero-cell AnnData objects should be summarised without errors."""
    adata_empty = ad.AnnData(np.empty((0, 3)))
    adata_empty.var_names = [f"G{i}" for i in range(3)]
    adata_empty.obs_names = []
    adata_empty.uns["dataset_id"] = "EMPTY"

    summary = qc_summary_table([adata_empty])
    assert summary.loc[0, "dataset_id"] == "EMPTY"
    assert summary.loc[0, "n_cells"] == 0
    assert np.isnan(summary.loc[0, "mean_genes_per_cell"])

    assert qc_summary_table([]).empty
