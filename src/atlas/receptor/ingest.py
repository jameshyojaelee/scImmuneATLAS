"""Ingest immune receptor repertoires and harmonize annotations."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import pandas as pd

from ..utils import ensure_dir
from .config import dataset_receptor_config

DEFAULT_PARQUET_DIR = Path("data/interim")
DEFAULT_METADATA_DIR = Path("processed/metrics/tcr")


def get_output_path(dataset_id: str, *, directory: Path = DEFAULT_PARQUET_DIR) -> Path:
    return directory / f"{dataset_id}_receptor.parquet"


def get_metadata_path(dataset_id: str, *, directory: Path = DEFAULT_METADATA_DIR) -> Path:
    return directory / f"{dataset_id}_ingest.json"


def _resolved_dirs(dataset_cfg: Dict, config: Dict) -> Tuple[Path, Path]:
    receptor_cfg = dataset_receptor_config(dataset_cfg, config)
    outdir = Path(receptor_cfg.get("interim_dir", DEFAULT_PARQUET_DIR))
    metrics_dir = Path(
        receptor_cfg.get(
            "qc_metrics_dir",
            receptor_cfg.get("metrics_dir", DEFAULT_METADATA_DIR),
        )
    )
    return outdir, metrics_dir


def resolved_output_path(dataset_cfg: Dict, config: Dict) -> Path:
    outdir, _ = _resolved_dirs(dataset_cfg, config)
    return get_output_path(dataset_cfg["id"], directory=outdir)


def resolved_metadata_path(dataset_cfg: Dict, config: Dict) -> Path:
    _, metrics_dir = _resolved_dirs(dataset_cfg, config)
    return get_metadata_path(dataset_cfg["id"], directory=metrics_dir)


def _normalise_productive(value: object) -> Optional[bool]:
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"true", "t", "yes", "y", "productive", "1"}:
        return True
    if text in {"false", "f", "no", "n", "unproductive", "0", "non-productive"}:
        return False
    return None


def _infer_locus(record: Dict[str, object]) -> Optional[str]:
    for key in ("locus", "chain", "chain_type"):
        value = record.get(key)
        if value is None:
            continue
        token = str(value).upper()
        if token.startswith("TR"):
            return token[:3]
        if token.startswith("IG"):
            return token[:3]
    return None


def _resolve_barcode(raw: object, dataset_id: str, receptor_cfg: Dict) -> str:
    barcode = str(raw).strip()
    for suffix in receptor_cfg.get("strip_suffixes", []):
        if barcode.endswith(suffix):
            barcode = barcode[: -len(suffix)]
    prefix = receptor_cfg.get("barcode_prefix", "")
    suffix = receptor_cfg.get("barcode_suffix", "")
    barcode = f"{prefix}{barcode}{suffix}"
    if receptor_cfg.get("append_dataset_id", False):
        separator = receptor_cfg.get("dataset_separator", "_")
        barcode = f"{dataset_id}{separator}{barcode}"
    if receptor_cfg.get("barcode_upper", False):
        barcode = barcode.upper()
    if receptor_cfg.get("barcode_lower", False):
        barcode = barcode.lower()
    return barcode


def _read_table(path: Path, receptor_cfg: Dict) -> pd.DataFrame:
    fmt = receptor_cfg.get("format", "10x_vdj").lower()
    sep = receptor_cfg.get("delimiter")
    if sep is None:
        if fmt in {"airr", "tsv"}:
            sep = "\t"
        else:
            sep = ","
    logging.info("Loading receptor file %s (format=%s, sep=%r)", path, fmt, sep)
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep=sep)
    if fmt == "airr":
        # AIRR format stores sequences per row with required columns
        pass
    return df


def _column(df: pd.DataFrame, names: Iterable[str]) -> Optional[pd.Series]:
    for name in names:
        if name in df.columns:
            return df[name]
    return None


def normalise_table(
    dataset_id: str,
    df: pd.DataFrame,
    receptor_cfg: Dict,
) -> pd.DataFrame:
    """Return harmonised dataframe with canonical receptor columns."""
    barcode_series = _column(df, [
        receptor_cfg.get("barcode_col", "barcode"),
        "barcode",
        "cell_id",
        "cell_id_num",
    ])
    if barcode_series is None:
        raise ValueError("Unable to locate barcode column in receptor table")

    clonotype_series = _column(df, [
        receptor_cfg.get("clonotype_col", "clonotype_id"),
        "clonotype",
        "raw_clonotype_id",
        "clonotype_id",
    ])

    cdr3_nt_series = _column(df, [
        receptor_cfg.get("cdr3_nt_col", "cdr3_nt"),
        "cdr3",
        "cdr3_nt",
    ])
    cdr3_aa_series = _column(df, [
        receptor_cfg.get("cdr3_aa_col", "cdr3_aa"),
        "cdr3_aa",
        "cdr3_aa_seq",
    ])

    chain_series = _column(df, [
        receptor_cfg.get("chain_col", "chain"),
        "locus",
        "gene" ,
    ])
    locus_series = _column(df, [
        receptor_cfg.get("locus_col", "locus"),
        "chain",
        "chain_type",
    ])

    v_gene_series = _column(df, [receptor_cfg.get("v_gene_col", "v_gene"), "v_gene", "v_call"])
    d_gene_series = _column(df, [receptor_cfg.get("d_gene_col", "d_gene"), "d_gene", "d_call"])
    j_gene_series = _column(df, [receptor_cfg.get("j_gene_col", "j_gene"), "j_gene", "j_call"])

    productive_series = _column(df, [receptor_cfg.get("productivity_col", "productive"), "productive"])
    umis_series = _column(df, [receptor_cfg.get("umis_col", "umis"), "umis", "umi_count"])
    reads_series = _column(df, [receptor_cfg.get("reads_col", "reads"), "reads", "read_count"])

    barcode_values = barcode_series.fillna("").astype(str).str.strip()
    harmonised = pd.DataFrame({
        "dataset_id": dataset_id,
        "barcode": barcode_values,
    })
    harmonised["cell_id"] = harmonised["barcode"].apply(
        lambda raw: _resolve_barcode(raw, dataset_id, receptor_cfg)
    )
    if clonotype_series is not None:
        harmonised["clonotype_id"] = clonotype_series.fillna("").astype(str)
    else:
        harmonised["clonotype_id"] = ""
    if chain_series is not None:
        harmonised["chain"] = chain_series.fillna("").astype(str).str.upper()
    else:
        harmonised["chain"] = (
            locus_series.fillna("").astype(str).str.upper()
            if locus_series is not None
            else ""
        )
    if locus_series is not None:
        harmonised["locus"] = locus_series.fillna("").astype(str).str.upper()
    else:
        harmonised["locus"] = harmonised.apply(lambda row: _infer_locus(row), axis=1)

    if cdr3_nt_series is not None:
        harmonised["cdr3_nt"] = cdr3_nt_series.fillna("").astype(str)
        harmonised["cdr3_nt_length"] = harmonised["cdr3_nt"].map(len)
    else:
        harmonised["cdr3_nt"] = ""
        harmonised["cdr3_nt_length"] = 0

    if cdr3_aa_series is not None:
        harmonised["cdr3_aa"] = cdr3_aa_series.fillna("").astype(str)
        harmonised["cdr3_aa_length"] = harmonised["cdr3_aa"].map(len)
    else:
        harmonised["cdr3_aa"] = ""
        harmonised["cdr3_aa_length"] = 0

    harmonised["v_gene"] = (
        v_gene_series.fillna("").astype(str) if v_gene_series is not None else ""
    )
    harmonised["d_gene"] = (
        d_gene_series.fillna("").astype(str) if d_gene_series is not None else ""
    )
    harmonised["j_gene"] = (
        j_gene_series.fillna("").astype(str) if j_gene_series is not None else ""
    )

    if productive_series is not None:
        productive_norm = productive_series.apply(_normalise_productive)
        harmonised["productive"] = productive_norm.fillna(False)
    else:
        harmonised["productive"] = pd.Series([None] * len(harmonised))

    if umis_series is not None:
        harmonised["umis"] = pd.to_numeric(umis_series, errors="coerce")
    else:
        harmonised["umis"] = None
    if reads_series is not None:
        harmonised["reads"] = pd.to_numeric(reads_series, errors="coerce")
    else:
        harmonised["reads"] = None

    harmonised["source_row"] = df.index

    harmonised = harmonised.drop_duplicates(subset=["cell_id", "chain", "cdr3_nt", "cdr3_aa"])
    return harmonised.reset_index(drop=True)


def ingest_dataset(dataset_cfg: Dict, config: Dict) -> Path:
    dataset_id = dataset_cfg["id"]
    receptor_cfg = dataset_receptor_config(dataset_cfg, config)
    source = receptor_cfg.get("contigs") or receptor_cfg.get("path")
    if not source:
        raise ValueError(f"Dataset '{dataset_id}' receptor config missing 'contigs' path")

    contig_path = Path(source)
    df_raw = _read_table(contig_path, receptor_cfg)
    df = normalise_table(dataset_id, df_raw, receptor_cfg)

    outdir, metrics_dir = _resolved_dirs(dataset_cfg, config)
    ensure_dir(outdir)
    out_path = get_output_path(dataset_id, directory=outdir)
    df.to_parquet(out_path, index=False)

    ensure_dir(metrics_dir)
    metadata_path = get_metadata_path(dataset_id, directory=metrics_dir)
    summary = {
        "dataset_id": dataset_id,
        "source": str(contig_path),
        "n_rows": int(len(df_raw)),
        "n_harmonised": int(len(df)),
        "n_cells": int(df["cell_id"].nunique()),
        "n_clonotypes": int(df["clonotype_id"].replace("", pd.NA).dropna().nunique()),
    }
    metadata_path.write_text(json.dumps(summary, indent=2))
    logging.info("Wrote harmonised receptor table to %s", out_path)
    return out_path
