"""Receptor repertoire module orchestration for scImmuneATLAS."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Dict, Iterable, Optional

from ..utils import ensure_dir, load_config, setup_logging
from . import analytics, ingest, qc, viz, report as report_mod
from .config import (
    allow_missing as cfg_allow_missing,
    dataset_has_receptor,
    global_receptor_config,
    is_enabled as cfg_is_enabled,
)


def _receptor_cfg(config: Dict) -> Dict:
    return global_receptor_config(config)


def is_enabled(config: Dict) -> bool:
    """Return True when the receptor module is enabled in the config."""
    return cfg_is_enabled(config)


def allow_missing(config: Dict) -> bool:
    return cfg_allow_missing(config)


def datasets_by_id(config: Dict) -> Dict[str, Dict]:
    return {entry["id"]: entry for entry in config.get("datasets", [])}


def datasets_with_receptor(config: Dict) -> Iterable[str]:
    for entry in config.get("datasets", []):
        if dataset_has_receptor(entry) or entry.get("receptor"):
            yield entry["id"]


def ingest_dataset(dataset_id: str, config: Dict) -> Optional[Path]:
    dataset_map = datasets_by_id(config)
    dataset_cfg = dataset_map.get(dataset_id)
    if not dataset_cfg:
        raise ValueError(f"Dataset '{dataset_id}' not defined in configuration")
    if not (dataset_has_receptor(dataset_cfg) or dataset_cfg.get("receptor")):
        if allow_missing(config):
            logging.info("Skipping receptor ingest for %s (no receptor config)", dataset_id)
            return None
        raise ValueError(
            f"Dataset '{dataset_id}' has no receptor configuration but receptor module is required"
        )
    return ingest.ingest_dataset(dataset_cfg, config)


def qc_dataset(dataset_id: str, config: Dict) -> Optional[Path]:
    dataset_map = datasets_by_id(config)
    dataset_cfg = dataset_map.get(dataset_id)
    if not dataset_cfg or not (dataset_has_receptor(dataset_cfg) or dataset_cfg.get("receptor")):
        if allow_missing(config):
            logging.info("Skipping receptor QC for %s (no receptor config)", dataset_id)
            return None
        raise ValueError(
            f"Dataset '{dataset_id}' unavailable or missing receptor config for QC stage"
        )
    return qc.run_dataset_qc(dataset_cfg, config)


def run_all(config: Dict) -> None:
    """Execute analytics, figures, and report integration for the receptor module."""
    if not is_enabled(config):
        logging.info("Receptor module disabled; skipping")
        return

    dataset_ids = list(datasets_with_receptor(config))
    if not dataset_ids:
        if allow_missing(config):
            logging.warning(
                "Receptor module enabled but no dataset has receptor configuration; nothing to do"
            )
            return
        raise ValueError("Receptor module enabled but no datasets provide receptor information")

    dataset_map = datasets_by_id(config)
    for dataset_id in dataset_ids:
        dataset_cfg = dataset_map[dataset_id]
        ingest_path = ingest.resolved_output_path(dataset_cfg, config)
        if not ingest_path.exists():
            logging.info("Re-ingesting receptor data for %s", dataset_id)
            ingest_dataset(dataset_id, config)
        qc_path = qc.get_output_path(dataset_id, config)
        if not qc_path.exists():
            logging.info("Recomputing receptor QC metrics for %s", dataset_id)
            qc_dataset(dataset_id, config)

    analytics_outputs = analytics.run_analytics(config)
    viz.generate_figures(config, analytics_outputs)
    report_mod.write_report_fragments(config, analytics_outputs)


def run_pipeline(config: Dict) -> None:
    """Full receptor pipeline covering ingest, QC, analytics, and viz."""
    if not is_enabled(config):
        logging.info("Receptor module disabled; skipping")
        return

    for dataset_id in datasets_with_receptor(config):
        ingest_dataset(dataset_id, config)
        qc_dataset(dataset_id, config)

    run_all(config)


def _cli_entry(args: argparse.Namespace) -> None:
    setup_logging(args.log_level)
    config = load_config(args.config)

    if args.stage in {None, "all"}:
        if is_enabled(config) or args.force:
            run_pipeline(config)
        else:
            logging.info("Receptor module disabled. Use --force to override.")
        return

    if args.stage == "ingest":
        if not args.dataset:
            raise SystemExit("--dataset is required for ingest stage")
        ingest_dataset(args.dataset, config)
        return

    if args.stage == "qc":
        if not args.dataset:
            raise SystemExit("--dataset is required for qc stage")
        qc_dataset(args.dataset, config)
        return

    if args.stage == "analytics":
        run_all(config)
        return

    if args.stage == "viz":
        outputs = analytics.run_analytics(config)
        viz.generate_figures(config, outputs)
        return

    raise SystemExit(f"Unknown receptor stage: {args.stage}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Receptor repertoire pipeline")
    parser.add_argument("--config", default="config/atlas.yaml", help="Config file path")
    parser.add_argument(
        "--stage",
        choices=["ingest", "qc", "analytics", "viz", "all"],
        default="all",
        help="Stage to execute",
    )
    parser.add_argument("--dataset", help="Dataset identifier (for ingest/qc)")
    parser.add_argument("--log-level", default="INFO", help="Logging level")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Run even when receptor module disabled in configuration",
    )
    args = parser.parse_args()
    _cli_entry(args)


if __name__ == "__main__":
    main()
