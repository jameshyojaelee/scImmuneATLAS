"""Unit tests for receptor I/O helpers."""

from __future__ import annotations

import tempfile
from pathlib import Path

import pandas as pd

from src.atlas.receptor.config import dataset_receptor_config
from src.atlas.receptors_io import load_receptor_table


def test_dataset_receptor_config_prefers_dataset_fields():
    config = {
        "tcr": {
            "enabled": True,
            "metrics_dir": "processed/metrics/tcr",
            "figures_dir": "processed/figures/tcr",
            "format": "airr",
        }
    }
    dataset_cfg = {
        "id": "DATASET_A",
        "receptor_path": "data/raw/custom.csv",
        "receptor_format": "10x_vdj",
        "receptor_sha256": "abc123",
        "receptor": {"barcode_prefix": "CELL_"},
    }

    merged = dataset_receptor_config(dataset_cfg, config)

    assert merged["contigs"] == "data/raw/custom.csv"
    assert merged["format"] == "10x_vdj"  # dataset override
    assert merged["contigs_sha256"] == "abc123"
    # fallback from global config still preserved
    assert merged["metrics_dir"] == "processed/metrics/tcr"
    assert merged["barcode_prefix"] == "CELL_"


def test_load_receptor_table_10x(tmp_path: Path):
    csv_path = tmp_path / "contigs.csv"
    pd.DataFrame(
        {
            "barcode": ["AAACCTG_1", "AAACCTG_1"],
            "chain": ["TRA", "TRB"],
            "clonotype_id": ["clon1", "clon1"],
            "cdr3_nt": ["ATGC", "ATGA"],
            "cdr3_aa": ["MK", "ML"],
            "v_gene": ["TRAV1", "TRBV1"],
            "productive": ["True", "False"],
            "umis": [5, 7],
        }
    ).to_csv(csv_path, index=False)

    dataset_info = {
        "id": "DATASET_A",
        "receptor_path": str(csv_path),
        "receptor_format": "10x_vdj",
    }

    df = load_receptor_table(dataset_info, config={"tcr": {"enabled": True}})

    assert df.shape[0] == 2
    assert set(df.columns).issuperset({
        "dataset_id",
        "cell_id",
        "chain",
        "clonotype_id",
        "cdr3_nt",
        "cdr3_aa",
        "productive",
    })
    assert df["dataset_id"].unique().tolist() == ["DATASET_A"]
    # barcodes should remain strings without whitespace
    assert all(df["cell_id"].str.contains("AAACCTG"))
    # productive coerced to boolean
    assert df["productive"].dtype.name == "boolean"


def test_load_receptor_table_airr(tmp_path: Path):
    tsv_path = tmp_path / "airr.tsv"
    pd.DataFrame(
        {
            "cell_id": ["CELL1", "CELL1"],
            "locus": ["TRB", "TRB"],
            "cdr3": ["ATGC", "ATGA"],
            "cdr3_aa": ["MK", "ML"],
            "v_call": ["TRBV1", "TRBV2"],
            "productive": ["T", "F"],
        }
    ).to_csv(tsv_path, sep="\t", index=False)

    dataset_info = {
        "id": "DATASET_B",
        "receptor_path": str(tsv_path),
        "receptor_format": "airr",
    }

    df = load_receptor_table(dataset_info, config={"tcr": {"enabled": True}})

    assert df.shape[0] == 2
    assert df["chain"].str.upper().eq("TRB").all()
    assert df["v_gene"].iloc[0] == "TRBV1"
    assert df["cdr3_nt"].iloc[0] == "ATGC"


def test_jaccard_index():
    """Test Jaccard index calculation."""
    from src.atlas.tcr import _calculate_jaccard_index

    set1 = {"A", "B", "C"}
    set2 = {"B", "C", "D"}

    # Intersection: {B, C} = 2, Union: {A, B, C, D} = 4
    assert _calculate_jaccard_index(set1, set2) == 0.5

    # Identical sets
    assert _calculate_jaccard_index(set1, set1) == 1.0

    # No overlap
    set3 = {"X", "Y", "Z"}
    assert _calculate_jaccard_index(set1, set3) == 0.0

    # Empty sets
    assert _calculate_jaccard_index(set(), set()) == 0.0


def test_morisita_horn():
    """Test Morisita-Horn similarity calculation."""
    from collections import Counter
    from src.atlas.tcr import _calculate_morisita_horn

    counts1 = Counter({"A": 10, "B": 5, "C": 3})
    counts2 = Counter({"A": 8, "B": 6, "C": 2})

    # Should return a value between 0 and 1
    similarity = _calculate_morisita_horn(counts1, counts2)
    assert 0.0 <= similarity <= 1.0

    # Identical distributions
    assert _calculate_morisita_horn(counts1, counts1) == 1.0

    # No overlap
    counts3 = Counter({"X": 5, "Y": 3})
    assert _calculate_morisita_horn(counts1, counts3) == 0.0

    # Empty counters
    assert _calculate_morisita_horn(Counter(), Counter()) == 0.0


def test_repertoire_overlap():
    """Test repertoire overlap computation."""
    from src.atlas.tcr import _compute_repertoire_overlap

    datasets_receptor = {
        "DATASET_A": pd.DataFrame({
            "clonotype_id": ["clone1", "clone2", "clone3", "clone3"],
        }),
        "DATASET_B": pd.DataFrame({
            "clonotype_id": ["clone2", "clone3", "clone4", "clone4"],
        }),
    }

    result = _compute_repertoire_overlap(datasets_receptor)

    assert "dataset_ids" in result
    assert "jaccard" in result
    assert "morisita_horn" in result
    assert len(result["dataset_ids"]) == 2
    assert len(result["jaccard"]) == 2
    assert len(result["jaccard"][0]) == 2

    # Diagonal should be 1.0 (self-similarity)
    assert result["jaccard"][0][0] == 1.0
    assert result["jaccard"][1][1] == 1.0
    assert result["morisita_horn"][0][0] == 1.0
    assert result["morisita_horn"][1][1] == 1.0


def test_public_clonotypes():
    """Test public clonotype detection."""
    from src.atlas.tcr import _identify_public_clonotypes

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

    result = _identify_public_clonotypes(datasets_receptor, min_datasets=2)

    assert "n_public_clonotypes" in result
    assert "top_public_clonotypes" in result
    assert "sharing_distribution" in result

    # clone2 is shared across all 3 datasets
    assert result["n_public_clonotypes"] >= 1

    # Check that clone2 is in top public clonotypes
    top_clones = [c["clonotype_id"] for c in result["top_public_clonotypes"]]
    assert "clone2" in top_clones

    # clone2 should be in 3 datasets
    clone2_info = next(c for c in result["top_public_clonotypes"] if c["clonotype_id"] == "clone2")
    assert clone2_info["n_datasets"] == 3
    assert set(clone2_info["datasets"]) == {"DATASET_A", "DATASET_B", "DATASET_C"}


def test_cdr3_properties():
    """Test CDR3 property analysis."""
    from src.atlas.tcr import _analyze_cdr3_properties

    receptor_df = pd.DataFrame({
        "cdr3_aa": ["CASSRG", "CASSLK", "CASSDE", "CASSRN", ""],
    })

    result = _analyze_cdr3_properties(receptor_df)

    assert "mean_length" in result
    assert "median_length" in result
    assert "mean_charge" in result
    assert "mean_hydrophobicity" in result
    assert "length_distribution" in result

    # Mean length should be 6 (excluding empty string)
    assert result["mean_length"] == 6.0
    assert result["median_length"] == 6.0

    # Should have charge and hydrophobicity values
    assert isinstance(result["mean_charge"], float)
    assert isinstance(result["mean_hydrophobicity"], float)


def test_cdr3_properties_empty():
    """Test CDR3 property analysis with empty data."""
    from src.atlas.tcr import _analyze_cdr3_properties

    receptor_df = pd.DataFrame({
        "cdr3_aa": ["", "", ""],
    })

    result = _analyze_cdr3_properties(receptor_df)

    assert result["mean_length"] == 0
    assert result["mean_charge"] == 0
    assert result["mean_hydrophobicity"] == 0
