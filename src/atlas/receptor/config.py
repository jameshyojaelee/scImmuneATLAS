"""Configuration helpers for receptor/TCR module."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

DEFAULT_METRICS_DIR = "processed/metrics/tcr"
DEFAULT_FIGURES_DIR = "processed/figures/tcr"
DEFAULT_ATLAS_PATH = "processed/integrated_annotated.h5ad"


def _base_config(config: Optional[Dict]) -> Dict:
    if not config:
        return {}
    if "tcr" in config and config.get("tcr") is not None:
        return dict(config.get("tcr") or {})
    if "receptor" in config and config.get("receptor") is not None:
        return dict(config.get("receptor") or {})
    return {}


def global_receptor_config(config: Optional[Dict]) -> Dict:
    cfg = _base_config(config)
    cfg.setdefault("qc_metrics_dir", cfg.get("metrics_dir", DEFAULT_METRICS_DIR))
    cfg.setdefault("figures_dir", cfg.get("figures_dir", DEFAULT_FIGURES_DIR))
    cfg.setdefault("atlas_path", cfg.get("atlas_path", DEFAULT_ATLAS_PATH))
    cfg.setdefault("allow_missing", bool(cfg.get("allow_missing", False)))
    cfg.setdefault("enabled", bool(cfg.get("enabled", False)))
    return cfg


def is_enabled(config: Optional[Dict]) -> bool:
    return bool(global_receptor_config(config).get("enabled", False))


def allow_missing(config: Optional[Dict]) -> bool:
    return bool(global_receptor_config(config).get("allow_missing", False))


def dataset_has_receptor(dataset_cfg: Dict) -> bool:
    if dataset_cfg.get("receptor_path"):
        return True
    nested = dataset_cfg.get("receptor") or {}
    return bool(nested.get("contigs") or nested.get("path") or nested.get("url"))


def dataset_receptor_config(dataset_cfg: Dict, config: Optional[Dict]) -> Dict:
    cfg = global_receptor_config(config)
    nested = dataset_cfg.get("receptor") or {}
    combined = cfg.copy()
    combined.update(nested)

    if dataset_cfg.get("receptor_path"):
        combined["contigs"] = dataset_cfg["receptor_path"]
    if dataset_cfg.get("receptor_format"):
        combined["format"] = dataset_cfg["receptor_format"]
    if dataset_cfg.get("receptor_sha256"):
        combined["contigs_sha256"] = dataset_cfg["receptor_sha256"]

    return combined


def metrics_dir(config: Optional[Dict]) -> Path:
    return Path(global_receptor_config(config).get("qc_metrics_dir", DEFAULT_METRICS_DIR))


def figures_dir(config: Optional[Dict]) -> Path:
    return Path(global_receptor_config(config).get("figures_dir", DEFAULT_FIGURES_DIR))
