"""Command-line interface for scImmuneATLAS."""

from __future__ import annotations

import argparse
import logging
from typing import Optional

from . import annotate, benchmark, doublets, integration, qc, export, viz
from .utils import generate_report, load_config, setup_logging


def _run_qc(config: dict, dataset: Optional[str]) -> None:
    if dataset:
        qc.process_dataset_qc(dataset, config)
    else:
        for dataset_info in config["datasets"]:
            qc.process_dataset_qc(dataset_info["id"], config)


def _run_doublets(config: dict, dataset: Optional[str]) -> None:
    if dataset:
        doublets.process_dataset_doublets(dataset, config)
    else:
        for dataset_info in config["datasets"]:
            doublets.process_dataset_doublets(dataset_info["id"], config)


def _run_integration(config: dict) -> None:
    integration.run_integration(config)


def _run_annotation(config: dict) -> None:
    annotate.run_annotation(config)


def _run_export(config: dict) -> None:
    export.run_export(config)


def _run_viz(config: dict) -> None:
    viz.generate_all_figures(config)


def _run_benchmark(config: dict) -> None:
    benchmark.run_benchmark(config)


def _run_report(config: dict, config_path: str) -> None:
    generate_report(config_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="scImmuneATLAS command-line interface")
    parser.add_argument("command", choices=[
        "qc",
        "doublets",
        "integrate",
        "annotate",
        "export",
        "viz",
        "benchmark",
        "report",
    ])
    parser.add_argument("--config", default="config/atlas.yaml", help="Config file path")
    parser.add_argument("--dataset", help="Dataset ID (for qc/doublets)")
    parser.add_argument("--log-level", default="INFO", help="Logging level")

    args = parser.parse_args()

    setup_logging(args.log_level)
    config = load_config(args.config)

    if args.command == "qc":
        _run_qc(config, args.dataset)
    elif args.command == "doublets":
        _run_doublets(config, args.dataset)
    elif args.command == "integrate":
        _run_integration(config)
    elif args.command == "annotate":
        _run_annotation(config)
    elif args.command == "export":
        _run_export(config)
    elif args.command == "viz":
        _run_viz(config)
    elif args.command == "benchmark":
        _run_benchmark(config)
    elif args.command == "report":
        _run_report(config, args.config)
    else:  # pragma: no cover - safeguard
        raise ValueError(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()
