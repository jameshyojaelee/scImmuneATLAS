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
