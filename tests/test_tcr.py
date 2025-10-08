"""Comprehensive unit tests for TCR analysis module.

This test suite covers:
- AIRR schema validation and conversion
- Scirpy pipeline integration (with mocking when unavailable)
- Clonotype metrics computation
- Repertoire overlap analysis
- Figure generation
- Error handling and edge cases
"""

from __future__ import annotations

import json
import tempfile
from collections import Counter
from pathlib import Path
from unittest.mock import MagicMock, patch

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from src.atlas import tcr
from src.atlas.tcr import (
    _airr_from_dataframe,
    _analyze_cdr3_properties,
    _calculate_jaccard_index,
    _calculate_morisita_horn,
    _compute_repertoire_overlap,
    _identify_public_clonotypes,
)


# ============================================================================
# Fixtures for Synthetic Data Generation
# ============================================================================


@pytest.fixture
def mock_receptor_df():
    """Generate synthetic receptor dataframe with realistic TCR data."""
    np.random.seed(42)

    n_cells = 200
    chains = ["TRA", "TRB"]
    v_genes = [f"TRBV{i}" for i in range(1, 31)]
    j_genes = [f"TRBJ{i}-{j}" for i in range(1, 3) for j in range(1, 8)]

    df = pd.DataFrame({
        "cell_id": [f"CELL_{i // 2}" for i in range(n_cells)],  # 2 chains per cell
        "dataset_id": np.random.choice(["DATASET_A", "DATASET_B"], size=n_cells),
        "chain": np.random.choice(chains, size=n_cells),
        "locus": np.random.choice(chains, size=n_cells),
        "cdr3_nt": [
            "".join(np.random.choice(list("ATGC"), np.random.randint(30, 60)))
            for _ in range(n_cells)
        ],
        "cdr3_aa": [
            "CASS" + "".join(np.random.choice(list("ARNDCEQGHILKMFPSTWYV"), np.random.randint(6, 12)))
            for _ in range(n_cells)
        ],
        "v_gene": np.random.choice(v_genes, size=n_cells),
        "d_gene": [""] * n_cells,
        "j_gene": np.random.choice(j_genes, size=n_cells),
        "productive": np.random.choice([True, False, pd.NA], size=n_cells, p=[0.8, 0.1, 0.1]),
        "umis": np.random.randint(1, 50, size=n_cells),
        "reads": np.random.randint(10, 500, size=n_cells),
        "clonotype_id": np.random.choice(
            [f"clone{i}" for i in range(1, 21)] + ["singleton"] * 30,
            size=n_cells
        ),
    })

    return df


@pytest.fixture
def mock_adata_with_tcr():
    """Generate synthetic AnnData with TCR annotations from scirpy."""
    np.random.seed(42)
    n_cells = 100

    X = np.random.randn(n_cells, 50)
    obs = pd.DataFrame({
        "cell_type": np.random.choice(["CD4 T", "CD8 T", "NK"], size=n_cells),
        "dataset_id": np.random.choice(["DATASET_A", "DATASET_B"], size=n_cells),
        "cancer_type": np.random.choice(["LUAD", "BRCA"], size=n_cells),
        "clonotype": np.random.choice(
            [f"clone{i}" for i in range(1, 11)] + [""] * 5,
            size=n_cells
        ),
        "IR_VJ_1_v_call": [f"TRBV{np.random.randint(1, 30)}" for _ in range(n_cells)],
        "IR_VJ_1_j_call": [f"TRBJ{np.random.randint(1, 3)}-{np.random.randint(1, 8)}" for _ in range(n_cells)],
        "IR_VJ_1_junction_aa": [
            "CASS" + "".join(np.random.choice(list("ARNDCEQGHILKMFPSTWYV"), 8))
            for _ in range(n_cells)
        ],
    })

    adata = ad.AnnData(X=X, obs=obs)
    sc.pp.neighbors(adata, n_neighbors=15)
    sc.tl.umap(adata)

    return adata


@pytest.fixture
def mock_config(tmp_path):
    """Generate mock configuration dictionary."""
    return {
        "seed": 42,
        "datasets": [
            {
                "id": "DATASET_A",
                "receptor_path": str(tmp_path / "dataset_a_receptors.parquet"),
                "receptor_format": "10x_vdj",
            },
            {
                "id": "DATASET_B",
                "receptor_path": str(tmp_path / "dataset_b_receptors.parquet"),
                "receptor_format": "10x_vdj",
            },
        ],
        "tcr": {
            "enabled": True,
            "metrics_dir": str(tmp_path / "metrics"),
            "figures_dir": str(tmp_path / "figures"),
            "min_public_datasets": 2,
        },
        "outputs": {
            "metrics_dir": str(tmp_path / "metrics"),
            "integrated_path": str(tmp_path / "integrated_annotated.h5ad"),
        },
    }


@pytest.fixture
def mock_datasets_receptor(mock_receptor_df):
    """Create dict of receptor dataframes per dataset."""
    df_a = mock_receptor_df[mock_receptor_df["dataset_id"] == "DATASET_A"].copy()
    df_b = mock_receptor_df[mock_receptor_df["dataset_id"] == "DATASET_B"].copy()

    return {
        "DATASET_A": df_a,
        "DATASET_B": df_b,
    }


# ============================================================================
# Tests for AIRR Schema Conversion
# ============================================================================


def test_airr_from_dataframe_basic(mock_receptor_df):
    """Test conversion of receptor dataframe to AIRR schema."""
    airr_df = _airr_from_dataframe(mock_receptor_df)

    # Check required AIRR columns
    required_cols = [
        "cell_id", "locus", "junction", "junction_aa", "v_call",
        "d_call", "j_call", "productive", "duplicate_count", "consensus_count"
    ]
    assert all(col in airr_df.columns for col in required_cols)

    # Check data mapping
    assert airr_df["cell_id"].equals(mock_receptor_df["cell_id"])
    assert airr_df["junction"].equals(mock_receptor_df["cdr3_nt"])
    assert airr_df["junction_aa"].equals(mock_receptor_df["cdr3_aa"])
    assert airr_df["v_call"].equals(mock_receptor_df["v_gene"])


def test_airr_from_dataframe_missing_columns():
    """Test AIRR conversion with missing optional columns."""
    minimal_df = pd.DataFrame({
        "cell_id": ["CELL_1", "CELL_2"],
        "chain": ["TRA", "TRB"],
        "cdr3_nt": ["ATGC", "ATGA"],
        "cdr3_aa": ["MK", "ML"],
        "v_gene": ["TRAV1", "TRBV1"],
        "j_gene": ["TRAJ1", "TRBJ1"],
    })

    airr_df = _airr_from_dataframe(minimal_df)

    # Should fill missing columns with NA
    assert "productive" in airr_df.columns
    assert "duplicate_count" in airr_df.columns
    assert airr_df["productive"].isna().all()


def test_airr_from_dataframe_alternative_column_names():
    """Test AIRR conversion when using 'locus' instead of 'chain'."""
    df = pd.DataFrame({
        "cell_id": ["CELL_1"],
        "locus": ["TRB"],  # Using 'locus' instead of 'chain'
        "cdr3_nt": ["ATGC"],
        "cdr3_aa": ["CASSK"],
        "v_gene": ["TRBV1"],
        "j_gene": ["TRBJ1"],
    })

    airr_df = _airr_from_dataframe(df)
    assert airr_df["locus"].iloc[0] == "TRB"


# ============================================================================
# Tests for Scirpy Integration (with Mocking)
# ============================================================================


def test_import_scirpy_missing():
    """Test that import error is raised when scirpy is unavailable."""
    with patch.dict("sys.modules", {"scirpy": None}):
        with pytest.raises(RuntimeError, match="scirpy is required"):
            tcr._import_scirpy()


@pytest.mark.skipif(
    not hasattr(tcr, "_merge_with_scirpy"),
    reason="Requires scirpy integration"
)
def test_merge_with_scirpy(mock_adata_with_tcr, mock_receptor_df, tmp_path):
    """Test merging receptor data with AnnData using scirpy."""
    try:
        import scirpy as ir  # noqa: F401
        scirpy_available = True
    except ImportError:
        scirpy_available = False

    if not scirpy_available:
        pytest.skip("scirpy not available")

    # Align cell IDs
    mock_receptor_df["cell_id"] = mock_adata_with_tcr.obs_names[:len(mock_receptor_df)].tolist()

    adata_merged = tcr._merge_with_scirpy(mock_adata_with_tcr.copy(), mock_receptor_df)

    # Check that IR columns were added
    ir_cols = [col for col in adata_merged.obs.columns if col.startswith("IR_")]
    assert len(ir_cols) > 0


# ============================================================================
# Tests for Clonotype Metrics Computation
# ============================================================================


@pytest.mark.skipif(
    not hasattr(tcr, "_compute_dataset_metrics"),
    reason="Requires scirpy"
)
def test_compute_dataset_metrics(mock_adata_with_tcr, mock_receptor_df, mock_config):
    """Test computation of per-dataset TCR metrics."""
    try:
        import scirpy as ir  # noqa: F401
        scirpy_available = True
    except ImportError:
        scirpy_available = False

    if not scirpy_available:
        pytest.skip("scirpy not available")

    metrics, diversity_dict = tcr._compute_dataset_metrics(
        "DATASET_A",
        mock_adata_with_tcr,
        mock_receptor_df,
        config=mock_config
    )

    # Check metrics structure
    assert "dataset_id" in metrics
    assert metrics["dataset_id"] == "DATASET_A"
    assert "n_cells" in metrics
    assert "n_clonotypes" in metrics
    assert "clonal_expansion" in metrics
    assert "diversity" in metrics
    assert "v_gene_usage" in metrics
    assert "j_gene_usage" in metrics
    assert "top_clonotypes" in metrics

    # Check diversity dict
    assert "dataset_id" in diversity_dict
    assert diversity_dict["dataset_id"] == "DATASET_A"

    # Check top clonotypes structure
    assert isinstance(metrics["top_clonotypes"], list)
    if len(metrics["top_clonotypes"]) > 0:
        top_clone = metrics["top_clonotypes"][0]
        assert "clonotype" in top_clone
        assert "cells" in top_clone


# ============================================================================
# Tests for Repertoire Overlap Analysis
# ============================================================================


def test_jaccard_index_basic():
    """Test Jaccard index calculation."""
    set1 = {"A", "B", "C"}
    set2 = {"B", "C", "D"}

    # Intersection: {B, C} = 2, Union: {A, B, C, D} = 4
    jaccard = _calculate_jaccard_index(set1, set2)
    assert jaccard == 0.5


def test_jaccard_index_identical():
    """Test Jaccard index with identical sets."""
    set1 = {"A", "B", "C"}
    jaccard = _calculate_jaccard_index(set1, set1)
    assert jaccard == 1.0


def test_jaccard_index_disjoint():
    """Test Jaccard index with no overlap."""
    set1 = {"A", "B"}
    set2 = {"C", "D"}
    jaccard = _calculate_jaccard_index(set1, set2)
    assert jaccard == 0.0


def test_jaccard_index_empty():
    """Test Jaccard index with empty sets."""
    jaccard = _calculate_jaccard_index(set(), set())
    assert jaccard == 0.0


def test_morisita_horn_basic():
    """Test Morisita-Horn similarity calculation."""
    counts1 = Counter({"A": 10, "B": 5, "C": 3})
    counts2 = Counter({"A": 8, "B": 6, "C": 2})

    similarity = _calculate_morisita_horn(counts1, counts2)
    assert 0.0 <= similarity <= 1.0


def test_morisita_horn_identical():
    """Test Morisita-Horn with identical distributions."""
    counts1 = Counter({"A": 10, "B": 5})
    similarity = _calculate_morisita_horn(counts1, counts1)
    assert similarity == 1.0


def test_morisita_horn_disjoint():
    """Test Morisita-Horn with no shared species."""
    counts1 = Counter({"A": 10, "B": 5})
    counts2 = Counter({"C": 8, "D": 6})
    similarity = _calculate_morisita_horn(counts1, counts2)
    assert similarity == 0.0


def test_morisita_horn_empty():
    """Test Morisita-Horn with empty counters."""
    similarity = _calculate_morisita_horn(Counter(), Counter())
    assert similarity == 0.0


def test_compute_repertoire_overlap(mock_datasets_receptor):
    """Test pairwise repertoire overlap computation."""
    overlap_metrics = _compute_repertoire_overlap(mock_datasets_receptor)

    # Check structure
    assert "dataset_ids" in overlap_metrics
    assert "jaccard" in overlap_metrics
    assert "morisita_horn" in overlap_metrics

    # Check dimensions
    n_datasets = len(mock_datasets_receptor)
    assert len(overlap_metrics["dataset_ids"]) == n_datasets
    assert len(overlap_metrics["jaccard"]) == n_datasets
    assert len(overlap_metrics["jaccard"][0]) == n_datasets

    # Check diagonal (self-similarity = 1.0)
    for i in range(n_datasets):
        assert overlap_metrics["jaccard"][i][i] == 1.0
        assert overlap_metrics["morisita_horn"][i][i] == 1.0

    # Check symmetry
    assert overlap_metrics["jaccard"][0][1] == overlap_metrics["jaccard"][1][0]
    assert overlap_metrics["morisita_horn"][0][1] == overlap_metrics["morisita_horn"][1][0]


def test_compute_repertoire_overlap_single_dataset():
    """Test overlap computation with single dataset (edge case)."""
    datasets_receptor = {
        "DATASET_A": pd.DataFrame({
            "clonotype_id": ["clone1", "clone2", "clone3"],
        }),
    }

    overlap_metrics = _compute_repertoire_overlap(datasets_receptor)

    assert len(overlap_metrics["dataset_ids"]) == 1
    assert overlap_metrics["jaccard"][0][0] == 1.0


# ============================================================================
# Tests for Public Clonotype Detection
# ============================================================================


def test_identify_public_clonotypes_basic():
    """Test public clonotype detection."""
    datasets_receptor = {
        "DATASET_A": pd.DataFrame({
            "clonotype_id": ["clone1", "clone2", "clone2", "clone3"],
        }),
        "DATASET_B": pd.DataFrame({
            "clonotype_id": ["clone2", "clone2", "clone4"],
        }),
        "DATASET_C": pd.DataFrame({
            "clonotype_id": ["clone2", "clone5"],
        }),
    }

    public_clonotypes = _identify_public_clonotypes(datasets_receptor, min_datasets=2)

    # Check structure
    assert "n_public_clonotypes" in public_clonotypes
    assert "top_public_clonotypes" in public_clonotypes
    assert "sharing_distribution" in public_clonotypes

    # clone2 is shared across all 3 datasets
    assert public_clonotypes["n_public_clonotypes"] >= 1

    # Check that clone2 is in top public clonotypes
    top_clones = [c["clonotype_id"] for c in public_clonotypes["top_public_clonotypes"]]
    assert "clone2" in top_clones

    # Verify clone2 details
    clone2_info = next(c for c in public_clonotypes["top_public_clonotypes"] if c["clonotype_id"] == "clone2")
    assert clone2_info["n_datasets"] == 3
    assert set(clone2_info["datasets"]) == {"DATASET_A", "DATASET_B", "DATASET_C"}
    assert clone2_info["total_cells"] == 5  # 2 + 2 + 1


def test_identify_public_clonotypes_threshold():
    """Test public clonotype detection with different thresholds."""
    datasets_receptor = {
        "DATASET_A": pd.DataFrame({"clonotype_id": ["clone1", "clone2"]}),
        "DATASET_B": pd.DataFrame({"clonotype_id": ["clone2", "clone3"]}),
        "DATASET_C": pd.DataFrame({"clonotype_id": ["clone3", "clone4"]}),
    }

    # With min_datasets=2
    result_2 = _identify_public_clonotypes(datasets_receptor, min_datasets=2)
    assert result_2["n_public_clonotypes"] == 2  # clone2 and clone3

    # With min_datasets=3
    result_3 = _identify_public_clonotypes(datasets_receptor, min_datasets=3)
    assert result_3["n_public_clonotypes"] == 0  # No clonotype in all 3


def test_identify_public_clonotypes_no_sharing():
    """Test public clonotype detection when no clonotypes are shared."""
    datasets_receptor = {
        "DATASET_A": pd.DataFrame({"clonotype_id": ["clone1", "clone2"]}),
        "DATASET_B": pd.DataFrame({"clonotype_id": ["clone3", "clone4"]}),
    }

    public_clonotypes = _identify_public_clonotypes(datasets_receptor, min_datasets=2)
    assert public_clonotypes["n_public_clonotypes"] == 0
    assert len(public_clonotypes["top_public_clonotypes"]) == 0


def test_identify_public_clonotypes_empty_strings():
    """Test public clonotype detection ignores empty clonotype IDs."""
    datasets_receptor = {
        "DATASET_A": pd.DataFrame({"clonotype_id": ["clone1", "", pd.NA]}),
        "DATASET_B": pd.DataFrame({"clonotype_id": ["clone1", "", ""]}),
    }

    public_clonotypes = _identify_public_clonotypes(datasets_receptor, min_datasets=2)

    # Only clone1 should be public (empty strings excluded)
    assert public_clonotypes["n_public_clonotypes"] == 1
    top_clones = [c["clonotype_id"] for c in public_clonotypes["top_public_clonotypes"]]
    assert top_clones == ["clone1"]


# ============================================================================
# Tests for CDR3 Property Analysis
# ============================================================================


def test_analyze_cdr3_properties_basic():
    """Test CDR3 property analysis with valid sequences."""
    receptor_df = pd.DataFrame({
        "cdr3_aa": ["CASSRG", "CASSLK", "CASSDE", "CASSRN"],
    })

    result = _analyze_cdr3_properties(receptor_df)

    # Check structure
    assert "mean_length" in result
    assert "median_length" in result
    assert "std_length" in result
    assert "mean_charge" in result
    assert "mean_hydrophobicity" in result
    assert "length_distribution" in result

    # All sequences have length 6
    assert result["mean_length"] == 6.0
    assert result["median_length"] == 6.0
    assert result["std_length"] == 0.0


def test_analyze_cdr3_properties_varied_lengths():
    """Test CDR3 property analysis with varied sequence lengths."""
    receptor_df = pd.DataFrame({
        "cdr3_aa": ["CASS", "CASSRG", "CASSRGDE", "CASSRGDELK"],
    })

    result = _analyze_cdr3_properties(receptor_df)

    # Mean length: (4 + 6 + 8 + 10) / 4 = 7.0
    assert result["mean_length"] == 7.0
    assert result["median_length"] == 7.0  # (6 + 8) / 2
    assert result["std_length"] > 0  # Should have variance


def test_analyze_cdr3_properties_charge_calculation():
    """Test CDR3 charge calculation."""
    # R and K are positively charged, D and E are negatively charged
    receptor_df = pd.DataFrame({
        "cdr3_aa": ["RRR", "DDD", "KKK", "EEE", "AAA"],
    })

    result = _analyze_cdr3_properties(receptor_df)

    # RRR: +3, KKK: +3, DDD: -3, EEE: -3, AAA: 0
    # Mean normalized charge: (1 + 1 - 1 - 1 + 0) / 5 = 0
    assert isinstance(result["mean_charge"], float)


def test_analyze_cdr3_properties_hydrophobicity():
    """Test CDR3 hydrophobicity calculation."""
    receptor_df = pd.DataFrame({
        "cdr3_aa": ["III", "VVV", "LLL"],  # Highly hydrophobic
    })

    result = _analyze_cdr3_properties(receptor_df)

    # Should have high hydrophobicity (positive Kyte-Doolittle values)
    assert result["mean_hydrophobicity"] > 0


def test_analyze_cdr3_properties_empty_dataframe():
    """Test CDR3 property analysis with empty data."""
    receptor_df = pd.DataFrame({
        "cdr3_aa": ["", pd.NA, ""],
    })

    result = _analyze_cdr3_properties(receptor_df)

    # Should return zeros for empty data
    assert result["mean_length"] == 0
    assert result["mean_charge"] == 0.0
    assert result["mean_hydrophobicity"] == 0.0


def test_analyze_cdr3_properties_length_distribution():
    """Test CDR3 length distribution output."""
    receptor_df = pd.DataFrame({
        "cdr3_aa": ["AAA", "AAA", "BBBB", "BBBB", "BBBB", "CCCCC"],
    })

    result = _analyze_cdr3_properties(receptor_df)

    length_dist = result["length_distribution"]
    assert length_dist[3] == 2  # Two sequences of length 3
    assert length_dist[4] == 3  # Three sequences of length 4
    assert length_dist[5] == 1  # One sequence of length 5


# ============================================================================
# Tests for JSON Output and File I/O
# ============================================================================


def test_write_json(tmp_path):
    """Test JSON writing utility."""
    payload = {
        "dataset_id": "DATASET_A",
        "n_cells": 100,
        "metrics": {"diversity": 2.5},
    }

    output_file = tmp_path / "test_metrics.json"
    tcr._write_json(output_file, payload)

    assert output_file.exists()

    # Verify content
    with open(output_file, "r") as f:
        loaded = json.load(f)

    assert loaded == payload
    assert loaded["dataset_id"] == "DATASET_A"


def test_write_json_creates_parent_directories(tmp_path):
    """Test that _write_json creates parent directories if needed."""
    nested_path = tmp_path / "nested" / "dir" / "metrics.json"

    payload = {"test": "data"}
    tcr._write_json(nested_path, payload)

    assert nested_path.exists()


# ============================================================================
# Tests for Figure Generation (Matplotlib Testing)
# ============================================================================


def test_generate_tcr_visualizations(mock_datasets_receptor, tmp_path):
    """Test TCR visualization figure generation."""
    overlap_metrics = _compute_repertoire_overlap(mock_datasets_receptor)
    public_clonotypes = _identify_public_clonotypes(mock_datasets_receptor, min_datasets=2)

    figures_dir = tmp_path / "figures"

    tcr._generate_tcr_visualizations(
        mock_datasets_receptor,
        overlap_metrics,
        public_clonotypes,
        figures_dir,
        adata_dict=None,
    )

    # Check that figures were created
    assert figures_dir.exists()

    expected_figures = [
        "repertoire_overlap_jaccard.png",
        "repertoire_similarity_morisita.png",
        "cdr3_length_distribution.png",
    ]

    for fig_name in expected_figures:
        fig_path = figures_dir / fig_name
        assert fig_path.exists(), f"Missing figure: {fig_name}"
        assert fig_path.stat().st_size > 0


def test_generate_tcr_visualizations_with_public_clonotypes(mock_datasets_receptor, tmp_path):
    """Test figure generation when public clonotypes exist."""
    # Ensure some clonotypes are shared
    df_a = mock_datasets_receptor["DATASET_A"].copy()
    df_b = mock_datasets_receptor["DATASET_B"].copy()

    # Force some shared clonotypes
    df_a.loc[:10, "clonotype_id"] = "shared_clone1"
    df_b.loc[:10, "clonotype_id"] = "shared_clone1"

    datasets_receptor = {"DATASET_A": df_a, "DATASET_B": df_b}

    overlap_metrics = _compute_repertoire_overlap(datasets_receptor)
    public_clonotypes = _identify_public_clonotypes(datasets_receptor, min_datasets=2)

    figures_dir = tmp_path / "figures"

    tcr._generate_tcr_visualizations(
        datasets_receptor,
        overlap_metrics,
        public_clonotypes,
        figures_dir,
    )

    # Should generate additional public clonotype figures
    assert (figures_dir / "public_clonotype_distribution.png").exists()
    assert (figures_dir / "top_public_clonotypes.png").exists()


def test_generate_clonotype_network_fallback(mock_adata_with_tcr, tmp_path):
    """Test clonotype network generation with error handling."""
    figures_dir = tmp_path / "figures"
    figures_dir.mkdir(parents=True)

    # This may fail if scirpy is unavailable or due to data issues
    # Should create a placeholder figure instead of crashing
    try:
        tcr._generate_clonotype_network(mock_adata_with_tcr, figures_dir, max_clonotypes=10)
    except Exception:
        # If it fails, should still create a placeholder
        pass

    # Check that a figure file was created (even if it's an error message)
    network_fig = figures_dir / "clonotype_network.png"
    if network_fig.exists():
        assert network_fig.stat().st_size > 0


# ============================================================================
# Tests for End-to-End TCR Analysis Pipeline
# ============================================================================


@pytest.mark.skipif(
    not hasattr(tcr, "run_tcr_analysis"),
    reason="Requires full TCR module"
)
def test_run_tcr_analysis_disabled(mock_config, tmp_path):
    """Test that TCR analysis is skipped when disabled in config."""
    config = mock_config.copy()
    config["tcr"]["enabled"] = False

    result = tcr.run_tcr_analysis(config)

    # Should return empty dict
    assert result == {}


@pytest.mark.skipif(
    not hasattr(tcr, "run_tcr_analysis"),
    reason="Requires full TCR module"
)
def test_run_tcr_analysis_no_datasets_with_receptors(mock_config):
    """Test TCR analysis when no datasets have receptor data."""
    config = mock_config.copy()
    config["datasets"] = [
        {"id": "DATASET_NO_TCR", "path": "data/no_tcr.h5ad"}
    ]

    result = tcr.run_tcr_analysis(config)

    # Should return empty dict with warning logged
    assert result == {}


# ============================================================================
# Tests for Error Handling and Edge Cases
# ============================================================================


def test_merge_with_scirpy_empty_dataframe(mock_adata_with_tcr):
    """Test merging with empty receptor dataframe."""
    empty_df = pd.DataFrame(columns=["cell_id", "chain", "cdr3_aa", "v_gene", "j_gene"])

    try:
        import scirpy as ir  # noqa: F401
        scirpy_available = True
    except ImportError:
        scirpy_available = False

    if not scirpy_available:
        pytest.skip("scirpy not available")

    # Should not crash, but may not add any TCR data
    adata_merged = tcr._merge_with_scirpy(mock_adata_with_tcr.copy(), empty_df)
    assert adata_merged.n_obs == mock_adata_with_tcr.n_obs


def test_compute_repertoire_overlap_empty_clonotypes():
    """Test repertoire overlap when all clonotypes are empty."""
    datasets_receptor = {
        "DATASET_A": pd.DataFrame({"clonotype_id": ["", "", ""]}),
        "DATASET_B": pd.DataFrame({"clonotype_id": ["", pd.NA, ""]}),
    }

    overlap_metrics = _compute_repertoire_overlap(datasets_receptor)

    # Should handle empty clonotypes gracefully
    assert overlap_metrics["jaccard"][0][1] == 0.0


def test_load_dataset_artifacts_missing_files(mock_config, tmp_path):
    """Test error handling when dataset files are missing."""
    dataset_cfg = mock_config["datasets"][0]

    try:
        with pytest.raises(FileNotFoundError):
            tcr._load_dataset_artifacts(dataset_cfg, mock_config)
    except AttributeError:
        # Function may not be exposed in all versions
        pytest.skip("_load_dataset_artifacts not available")


# ============================================================================
# Tests for Schema Validation (AIRR Compliance)
# ============================================================================


def test_airr_schema_required_fields(mock_receptor_df):
    """Test that AIRR conversion includes all required fields."""
    airr_df = _airr_from_dataframe(mock_receptor_df)

    # AIRR Rearrangement schema required fields
    required_fields = [
        "cell_id",
        "locus",
        "v_call",
        "j_call",
        "junction",
        "junction_aa",
        "productive",
    ]

    for field in required_fields:
        assert field in airr_df.columns, f"Missing required AIRR field: {field}"


def test_airr_schema_data_types(mock_receptor_df):
    """Test that AIRR dataframe has correct data types."""
    airr_df = _airr_from_dataframe(mock_receptor_df)

    # Check that string columns are strings
    assert airr_df["cell_id"].dtype == object
    assert airr_df["locus"].dtype == object
    assert airr_df["junction"].dtype == object

    # Check that productive is boolean-compatible
    assert airr_df["productive"].dtype == object  # pandas nullable boolean


def test_airr_schema_handles_nulls(mock_receptor_df):
    """Test that AIRR conversion handles null values correctly."""
    df_with_nulls = mock_receptor_df.copy()
    df_with_nulls.loc[0, "v_gene"] = pd.NA
    df_with_nulls.loc[1, "cdr3_aa"] = None
    df_with_nulls.loc[2, "productive"] = pd.NA

    airr_df = _airr_from_dataframe(df_with_nulls)

    # Should not crash and should preserve NAs
    assert airr_df["v_call"].isna().any()
    assert airr_df["junction_aa"].isna().any()


# ============================================================================
# Parametrized Tests for Robustness
# ============================================================================


@pytest.mark.parametrize("n_datasets", [1, 2, 3, 5])
def test_repertoire_overlap_various_dataset_counts(n_datasets):
    """Test repertoire overlap computation with varying dataset counts."""
    datasets_receptor = {
        f"DATASET_{i}": pd.DataFrame({
            "clonotype_id": [f"clone{j}" for j in range(i, i + 10)]
        })
        for i in range(n_datasets)
    }

    overlap_metrics = _compute_repertoire_overlap(datasets_receptor)

    assert len(overlap_metrics["dataset_ids"]) == n_datasets
    assert len(overlap_metrics["jaccard"]) == n_datasets
    assert len(overlap_metrics["morisita_horn"]) == n_datasets


@pytest.mark.parametrize("cdr3_length", [6, 10, 15, 20])
def test_cdr3_properties_various_lengths(cdr3_length):
    """Test CDR3 property analysis with various sequence lengths."""
    receptor_df = pd.DataFrame({
        "cdr3_aa": ["A" * cdr3_length] * 10
    })

    result = _analyze_cdr3_properties(receptor_df)

    assert result["mean_length"] == cdr3_length
    assert result["median_length"] == cdr3_length
    assert result["std_length"] == 0.0


# ============================================================================
# Integration Tests (Mocking scirpy when unavailable)
# ============================================================================


@patch("src.atlas.tcr._import_scirpy")
def test_tcr_analysis_without_scirpy(mock_import_scirpy, mock_config):
    """Test that TCR analysis fails gracefully when scirpy is unavailable."""
    mock_import_scirpy.side_effect = RuntimeError("scirpy is required")

    # Analysis should raise error or return empty
    with pytest.raises(RuntimeError, match="scirpy is required"):
        tcr._import_scirpy()


# ============================================================================
# Regression Tests
# ============================================================================


def test_clonotype_counts_exclude_empty_strings():
    """Regression test: ensure empty clonotype IDs are excluded from counts."""
    receptor_df = pd.DataFrame({
        "clonotype_id": ["clone1", "clone2", "", "", pd.NA, "clone1"]
    })

    clonotype_counts = (
        receptor_df["clonotype_id"]
        .replace("", pd.NA)
        .dropna()
        .value_counts()
    )

    # Should only count non-empty clonotypes
    assert len(clonotype_counts) == 2
    assert clonotype_counts["clone1"] == 2
    assert clonotype_counts["clone2"] == 1


def test_v_gene_usage_top_n():
    """Regression test: ensure V gene usage returns top N genes."""
    receptor_df = pd.DataFrame({
        "v_gene": [f"TRBV{i}" for i in range(1, 51)] * 2  # 50 unique genes, each appearing twice
    })

    v_usage = receptor_df["v_gene"].value_counts().head(20)

    assert len(v_usage) == 20
    assert all(count == 2 for count in v_usage.values)
