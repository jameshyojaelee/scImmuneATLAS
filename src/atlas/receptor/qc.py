"""Compute QA/QC metrics for immune receptor repertoires."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, Tuple

import pandas as pd

from ..utils import ensure_dir
from . import ingest
from .config import dataset_receptor_config, metrics_dir as config_metrics_dir

DEFAULT_QC_DIR = Path("processed/metrics/tcr")


def get_output_dir(config: Dict) -> Path:
    return config_metrics_dir(config)


def get_output_path(dataset_id: str, config: Dict) -> Path:
    return get_output_dir(config) / f"{dataset_id}_qc.json"


def _detect_modality(df: pd.DataFrame, dataset_cfg: Dict, config: Dict) -> str:
    receptor_cfg = dataset_receptor_config(dataset_cfg, config)
    modality = receptor_cfg.get("modality")
    if modality:
        return modality.lower()
    loci = df.get("locus")
    if loci is None:
        return "unknown"
    loci_set = set(str(val).upper() for val in loci.dropna())
    has_tr = any(loc.startswith("TR") for loc in loci_set)
    has_ig = any(loc.startswith("IG") for loc in loci_set)
    if has_tr and has_ig:
        return "both"
    if has_tr:
        return "tcr"
    if has_ig:
        return "bcr"
    return "unknown"


def _paired_label(chain_loci: set[str]) -> str:
    loci = {loc for loc in chain_loci if loc}
    if {"TRA", "TRB"}.issubset(loci):
        return "TCR_alpha_beta"
    if {"TRG", "TRD"}.issubset(loci):
        return "TCR_gamma_delta"
    if "IGH" in loci and ("IGK" in loci or "IGL" in loci):
        return "BCR_heavy_light"
    return "unpaired"


def _summarise_cells(df: pd.DataFrame) -> pd.DataFrame:
    grouped = df.groupby("cell_id")
    rows = []
    for cell_id, subset in grouped:
        loci = set(subset.get("locus", pd.Series(dtype=str)).dropna().astype(str))
        productive_flags = subset.get("productive", pd.Series(dtype=bool)).fillna(False)
        clonotypes = subset.get("clonotype_id", pd.Series(dtype=str)).replace("", pd.NA).dropna().unique()
        rows.append({
            "dataset_id": subset["dataset_id"].iloc[0],
            "cell_id": cell_id,
            "n_chains": int(len(subset)),
            "n_productive_chains": int(productive_flags.sum()),
            "has_productive_chain": bool(productive_flags.any()),
            "paired_label": _paired_label(loci),
            "clonotype_id": clonotypes[0] if len(clonotypes) else "",
        })
    return pd.DataFrame(rows)


def _aggregate_metrics(df: pd.DataFrame, cell_summary: pd.DataFrame, modality: str) -> Dict:
    productive_series = df.get("productive")
    productive_count = int(productive_series.fillna(False).sum()) if productive_series is not None else 0
    total_sequences = int(len(df))
    unique_cells = int(cell_summary["cell_id"].nunique()) if not cell_summary.empty else 0
    productive_cells = int(cell_summary["has_productive_chain"].sum()) if not cell_summary.empty else 0
    clonotypes = df.get("clonotype_id", pd.Series(dtype=str)).replace("", pd.NA).dropna()
    pairing_counts = cell_summary["paired_label"].value_counts().to_dict() if not cell_summary.empty else {}

    metrics = {
        "modality": modality,
        "n_sequences": total_sequences,
        "n_cells_with_receptor": unique_cells,
        "n_productive_sequences": productive_count,
        "pct_productive_sequences": (productive_count / total_sequences) if total_sequences else 0.0,
        "n_productive_cells": productive_cells,
        "pct_productive_cells": (productive_cells / unique_cells) if unique_cells else 0.0,
        "unique_clonotypes": int(clonotypes.nunique()),
        "median_chains_per_cell": float(cell_summary["n_chains"].median()) if not cell_summary.empty else 0.0,
        "mean_chains_per_cell": float(cell_summary["n_chains"].mean()) if not cell_summary.empty else 0.0,
        "paired_counts": pairing_counts,
    }
    return metrics


def run_dataset_qc(dataset_cfg: Dict, config: Dict) -> Path:
    dataset_id = dataset_cfg["id"]
    parquet_path = ingest.resolved_output_path(dataset_cfg, config)
    if not parquet_path.exists():
        raise FileNotFoundError(f"Receptor table missing for dataset {dataset_id}: {parquet_path}")

    df = pd.read_parquet(parquet_path)
    if df.empty:
        logging.warning("Receptor table for %s is empty", dataset_id)

    cell_summary = _summarise_cells(df) if not df.empty else pd.DataFrame()
    modality = _detect_modality(df, dataset_cfg, config)
    metrics = _aggregate_metrics(df, cell_summary, modality)

    metrics_dir = get_output_dir(config)
    ensure_dir(metrics_dir)
    qc_path = get_output_path(dataset_id, config)
    qc_path.write_text(json.dumps(metrics, indent=2))

    if not cell_summary.empty:
        cells_path = metrics_dir / f"{dataset_id}_cell_metrics.tsv"
        cell_summary.to_csv(cells_path, sep="\t", index=False)
    logging.info("Receptor QC metrics written to %s", qc_path)
    return qc_path
