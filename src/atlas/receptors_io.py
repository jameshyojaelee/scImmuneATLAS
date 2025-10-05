"""Helpers for loading and attaching immune receptor metadata during I/O."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import anndata as ad
import pandas as pd

from .schemas import validate_receptor_table

DEFAULT_CHAIN_TYPES = {
    "alpha_beta": {"TRA", "TRB"},
    "gamma_delta": {"TRG", "TRD"},
}

SUPPORTED_FORMATS = {"10x_vdj", "airr"}


def _tcr_config(config: Optional[Dict]) -> Dict:
    if not config:
        return {}
    if "tcr" in config:
        return config["tcr"] or {}
    if "receptor" in config:  # backwards compatibility
        return config["receptor"] or {}
    return {}


def _dataset_receptor_fields(dataset_info: Dict) -> Tuple[Optional[str], Optional[str]]:
    if not dataset_info:
        return None, None
    nested = dataset_info.get("receptor", {})
    path = (
        dataset_info.get("receptor_path")
        or nested.get("contigs")
        or nested.get("path")
    )
    fmt = dataset_info.get("receptor_format") or nested.get("format")
    return path, fmt


def _resolve_chain_types(config: Optional[Dict]) -> Dict[str, set[str]]:
    cfg = _tcr_config(config)
    chain_types = cfg.get("chain_types")
    resolved: Dict[str, set[str]] = {}
    if isinstance(chain_types, dict):
        resolved = {name: set(values or []) for name, values in chain_types.items()}
    elif isinstance(chain_types, list):
        for item in chain_types:
            if isinstance(item, dict) and "name" in item:
                resolved[item["name"]] = set(item.get("chains", []))
            elif isinstance(item, (list, tuple)) and len(item) >= 2:
                label = str(item[0])
                resolved[label] = set(item[1:])
    if not resolved:
        resolved = DEFAULT_CHAIN_TYPES.copy()
    return {name: {chain.upper() for chain in chains if chain} for name, chains in resolved.items()}


def _normalise_productive(series: pd.Series) -> pd.Series:
    mapping = {
        "true": True,
        "t": True,
        "yes": True,
        "y": True,
        "1": True,
        "false": False,
        "f": False,
        "no": False,
        "n": False,
        "0": False,
        "productive": True,
        "unproductive": False,
    }
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .map(mapping)
        .where(lambda x: x.isin([True, False]), other=pd.NA)
    )


def _load_10x_vdj(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.rename(
        columns={
            "barcode": "cell_id",
            "chain": "chain",
            "clonotype_id": "clonotype_id",
            "raw_clonotype_id": "clonotype_id",
            "cdr3_nt": "cdr3_nt",
            "cdr3": "cdr3_nt",
            "cdr3_aa": "cdr3_aa",
            "umis": "umis",
            "umi_count": "umis",
            "reads": "reads",
            "read_count": "reads",
            "v_gene": "v_gene",
            "d_gene": "d_gene",
            "j_gene": "j_gene",
            "productive": "productive",
        }
    )
    if "cell_id" not in df.columns and "barcode" in df.columns:
        df["cell_id"] = df["barcode"].astype(str)
    return df


def _load_airr(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    return df.rename(
        columns={
            "cell_id": "cell_id",
            "cell_id_full": "cell_id",
            "locus": "locus",
            "chain": "chain",
            "clone_id": "clonotype_id",
            "clonotype_id": "clonotype_id",
            "productive": "productive",
            "sequence_id": "sequence_id",
            "cdr3": "cdr3_nt",
            "cdr3_nt": "cdr3_nt",
            "cdr3_aa": "cdr3_aa",
            "v_call": "v_gene",
            "d_call": "d_gene",
            "j_call": "j_gene",
        }
    )


def load_receptor_table(
    dataset_info: Dict,
    *,
    config: Optional[Dict] = None,
) -> Optional[pd.DataFrame]:
    path_str, fmt = _dataset_receptor_fields(dataset_info)
    if not path_str:
        return None

    path = Path(path_str)
    if fmt is None:
        fmt = "10x_vdj" if path.suffix.lower() == ".csv" else "airr"
    fmt = fmt.lower()
    if fmt not in SUPPORTED_FORMATS:
        logging.warning("Unsupported receptor format '%s'; skipping", fmt)
        return None

    if fmt == "10x_vdj":
        df = _load_10x_vdj(path)
    elif fmt == "airr":
        df = _load_airr(path)
    else:  # pragma: no cover
        return None

    df = df.copy()
    dataset_id = dataset_info.get("id", "unknown_dataset")
    df["dataset_id"] = dataset_id

    if "cell_id" not in df.columns:
        raise ValueError(f"Receptor table for {dataset_id} missing 'cell_id' column")

    df["chain"] = df.get("chain", df.get("locus", "")).astype(str).str.upper()
    if "locus" not in df.columns:
        df["locus"] = df["chain"]

    if "productive" in df.columns:
        df["productive"] = _normalise_productive(df["productive"]).astype("boolean")
    else:
        df["productive"] = pd.Series([pd.NA] * len(df), dtype="boolean")

    required_columns = ["cdr3_nt", "cdr3_aa", "v_gene", "d_gene", "j_gene", "clonotype_id", "umis", "reads"]
    for column in required_columns:
        if column not in df.columns:
            df[column] = pd.NA

    df["cell_id"] = df["cell_id"].astype(str)
    return df


def _pairing_label(chains: Iterable[str], chain_map: Dict[str, set[str]]) -> str:
    chain_set = {chain for chain in chains if chain}
    for label, required in chain_map.items():
        if required.issubset(chain_set):
            return label
    return "incomplete"


def _summarise_receptor(df: pd.DataFrame, chain_map: Dict[str, set[str]]) -> pd.DataFrame:
    grouped = df.groupby("cell_id")
    summary = grouped.agg(
        n_chains=("cell_id", "count"),
        n_productive=("productive", lambda s: int(s.fillna(False).sum())),
        clonotype_id=("clonotype_id", lambda s: next((v for v in s if pd.notna(v) and v != ""), "")),
    )
    pairing = grouped["chain"].agg(lambda s: _pairing_label(s, chain_map))
    summary["pairing"] = pairing
    summary["has_productive_chain"] = summary["n_productive"] > 0
    return summary


def attach_receptor_data(
    adata: ad.AnnData,
    dataset_info: Dict,
    *,
    config: Optional[Dict] = None,
) -> None:
    """Attach receptor summaries to AnnData in-place if data are available."""

    tcr_config = _tcr_config(config)
    if not tcr_config.get("enabled", False):
        return

    receptor_df = load_receptor_table(dataset_info, config=config)
    if receptor_df is None or receptor_df.empty:
        return

    chain_map = _resolve_chain_types(config)
    chain_sets = [list(chains) for chains in chain_map.values() if chains]
    validate_receptor_table(receptor_df, chain_sets=chain_sets)

    summary = _summarise_receptor(receptor_df, chain_map)
    obs = adata.obs
    obs["tcr_chain_count"] = (
        summary["n_chains"].reindex(obs.index).fillna(0).astype(int)
    )
    obs["tcr_productive_chains"] = (
        summary["n_productive"].reindex(obs.index).fillna(0).astype(int)
    )
    obs["tcr_has_productive_chain"] = (
        summary["has_productive_chain"].reindex(obs.index).fillna(False)
    )
    obs["tcr_clonotype"] = (
        summary["clonotype_id"].reindex(obs.index).fillna("")
    )
    obs["tcr_pairing"] = (
        summary["pairing"].reindex(obs.index).fillna("absent")
    )

    adata.uns.setdefault("tcr", {})
    adata.uns["tcr"]["chain_types"] = {key: sorted(value) for key, value in chain_map.items()}
    adata.uns["tcr"].setdefault("tables", {})
    adata.uns["tcr"]["tables"][dataset_info.get("id", "dataset")] = receptor_df

    scirpy_settings = tcr_config.get("scirpy_settings") or {}
    if scirpy_settings:
        adata.uns["tcr"]["scirpy_settings"] = scirpy_settings

    logging.info(
        "Attached receptor metadata for %s (%d sequences, %d cells)",
        dataset_info.get("id"),
        len(receptor_df),
        summary.shape[0],
    )
