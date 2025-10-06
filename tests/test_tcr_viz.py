"""Unit tests for TCR visualization functions."""

from __future__ import annotations

import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from src.atlas import tcr_viz


@pytest.fixture
def mock_adata():
    """Create mock AnnData with clonotype annotations and UMAP."""
    np.random.seed(42)
    n_cells = 500

    # Create synthetic data
    X = np.random.randn(n_cells, 100)
    obs = pd.DataFrame({
        "clonotype": np.random.choice(
            ["clone1", "clone2", "clone3", "clone4", "clone5", "singleton"] * 10 + [""] * 50,
            size=n_cells
        ),
        "dataset_id": np.random.choice(["DATASET_A", "DATASET_B"], size=n_cells),
        "cancer_type": np.random.choice(["LUAD", "BRCA", "COAD"], size=n_cells),
    })

    adata = ad.AnnData(X=X, obs=obs)

    # Generate UMAP
    sc.pp.neighbors(adata, n_neighbors=15)
    sc.tl.umap(adata)

    return adata


@pytest.fixture
def mock_receptor_df():
    """Create mock receptor dataframe."""
    np.random.seed(42)

    chains = ["TRA", "TRB", "TRG", "TRD"]
    n_rows = 800

    df = pd.DataFrame({
        "cell_id": [f"CELL_{i}" for i in range(n_rows)],
        "dataset_id": np.random.choice(["DATASET_A", "DATASET_B"], size=n_rows),
        "chain": np.random.choice(chains, size=n_rows),
        "cdr3_nt": ["ATGC" * np.random.randint(10, 20) for _ in range(n_rows)],
        "cdr3_aa": ["CASS" + "".join(np.random.choice(list("ARNDCEQGHILKMFPSTWYV"), 6)) for _ in range(n_rows)],
        "v_gene": [f"TRBV{np.random.randint(1, 30)}" for _ in range(n_rows)],
        "j_gene": [f"TRBJ{np.random.randint(1, 3)}" for _ in range(n_rows)],
        "clonotype_id": np.random.choice(
            ["clone1", "clone2", "clone3", "clone4", "clone5"] * 20 + ["singleton"] * 100,
            size=n_rows
        ),
    })

    return df


@pytest.fixture
def mock_diversity_df():
    """Create mock diversity metrics dataframe."""
    return pd.DataFrame({
        "dataset_id": ["DATASET_A", "DATASET_B", "DATASET_C"],
        "cancer_type": ["LUAD", "BRCA", "LUAD"],
        "shannon_entropy": [2.5, 3.1, 2.8],
        "normalized_shannon_entropy": [0.7, 0.85, 0.75],
        "D50": [10, 15, 12],
    })


def test_plot_clonotype_frequency(mock_adata, tmp_path):
    """Test clonotype frequency barplot generation."""
    outpath = tmp_path / "clonotype_freq.png"

    tcr_viz.plot_clonotype_frequency(
        mock_adata,
        clonotype_col="clonotype",
        top_n=10,
        outpath=outpath,
        dataset_col="dataset_id",
    )

    assert outpath.exists()
    assert outpath.stat().st_size > 0


def test_plot_repertoire_diversity(mock_diversity_df, tmp_path):
    """Test repertoire diversity boxplot generation."""
    outpath = tmp_path / "diversity.png"

    tcr_viz.plot_repertoire_diversity(
        mock_diversity_df,
        groupby="cancer_type",
        outpath=outpath,
    )

    assert outpath.exists()
    assert outpath.stat().st_size > 0


def test_plot_umap_clonotype_size(mock_adata, tmp_path):
    """Test UMAP colored by clonotype size."""
    outpath = tmp_path / "umap_clonotype.png"

    tcr_viz.plot_umap_clonotype_size(
        mock_adata,
        clonotype_col="clonotype",
        umap_key="X_umap",
        outpath=outpath,
    )

    assert outpath.exists()
    assert outpath.stat().st_size > 0


def test_plot_spectratype(mock_receptor_df, tmp_path):
    """Test spectratype CDR3 length distribution plot."""
    outpath = tmp_path / "spectratype.png"

    tcr_viz.plot_spectratype(
        mock_receptor_df,
        cdr3_col="cdr3_nt",
        chain_col="chain",
        dataset_col="dataset_id",
        outpath=outpath,
        separate_chains=True,
    )

    assert outpath.exists()
    assert outpath.stat().st_size > 0


def test_plot_vj_pairing_heatmap(mock_receptor_df, tmp_path):
    """Test V-J pairing heatmap generation."""
    outpath = tmp_path / "vj_pairing.png"

    tcr_viz.plot_vj_pairing_heatmap(
        mock_receptor_df,
        v_col="v_gene",
        j_col="j_gene",
        top_n=10,
        outpath=outpath,
    )

    assert outpath.exists()
    assert outpath.stat().st_size > 0


def test_plot_clonal_expansion_summary(mock_adata, tmp_path):
    """Test clonal expansion summary plot."""
    outpath = tmp_path / "expansion_summary.png"

    tcr_viz.plot_clonal_expansion_summary(
        mock_adata,
        clonotype_col="clonotype",
        groupby="cancer_type",
        outpath=outpath,
    )

    assert outpath.exists()
    assert outpath.stat().st_size > 0


def test_generate_all_tcr_figures(mock_adata, mock_receptor_df, mock_diversity_df, tmp_path):
    """Test complete visualization suite generation."""
    figures_dir = tmp_path / "tcr_figures"

    tcr_viz.generate_all_tcr_figures(
        mock_adata,
        mock_receptor_df,
        mock_diversity_df,
        figures_dir,
        groupby="cancer_type",
    )

    # Check that multiple figures were created
    assert figures_dir.exists()
    figures = list(figures_dir.glob("*.png"))
    assert len(figures) >= 4  # At least 4 different figure types


def test_plot_clonotype_frequency_no_clonotypes(mock_adata, tmp_path):
    """Test graceful handling when no clonotypes present."""
    # Create adata without clonotypes
    adata_no_clonotypes = mock_adata.copy()
    adata_no_clonotypes.obs["clonotype"] = ""

    outpath = tmp_path / "no_clonotypes.png"

    # Should not raise error
    tcr_viz.plot_clonotype_frequency(
        adata_no_clonotypes,
        clonotype_col="clonotype",
        outpath=outpath,
    )


def test_plot_spectratype_no_chains(mock_receptor_df, tmp_path):
    """Test spectratype plot without chain separation."""
    outpath = tmp_path / "spectratype_no_chains.png"

    # Remove chain column
    df_no_chains = mock_receptor_df.drop(columns=["chain"])

    tcr_viz.plot_spectratype(
        df_no_chains,
        cdr3_col="cdr3_nt",
        chain_col="chain",
        outpath=outpath,
        separate_chains=False,
    )

    assert outpath.exists()


def test_plot_umap_downsampling(tmp_path):
    """Test UMAP plot with downsampling for large datasets."""
    np.random.seed(42)

    # Create large dataset
    n_cells = 100000
    X = np.random.randn(n_cells, 50)
    obs = pd.DataFrame({
        "clonotype": np.random.choice(["clone1", "clone2", "clone3"], size=n_cells),
    })

    adata_large = ad.AnnData(X=X, obs=obs)
    sc.pp.neighbors(adata_large, n_neighbors=15)
    sc.tl.umap(adata_large)

    outpath = tmp_path / "umap_downsampled.png"

    tcr_viz.plot_umap_clonotype_size(
        adata_large,
        clonotype_col="clonotype",
        outpath=outpath,
        max_points=10000,  # Force downsampling
    )

    assert outpath.exists()
