"""Command-line interface for scImmuneATLAS."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

from . import annotate, benchmark, doublets, export, integration, qc, receptor, tcr, viz, workflow
from .io import download_if_needed, load_matrix
from .schemas import validate_anndata
from .utils import generate_report, load_config, setup_logging, ensure_dir


def _run_qc(config: dict, dataset: Optional[str]) -> None:
    if dataset:
        qc.process_dataset_qc(dataset, config)
    else:
        for dataset_id in workflow.dataset_ids(config):
            qc.process_dataset_qc(dataset_id, config)


def _run_doublets(config: dict, dataset: Optional[str]) -> None:
    if dataset:
        doublets.process_dataset_doublets(dataset, config)
    else:
        for dataset_id in workflow.dataset_ids(config):
            doublets.process_dataset_doublets(dataset_id, config)


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


def _run_pipeline_cli(
    config: dict,
    config_path: str,
    jobs: Optional[int],
    force_python: bool,
) -> None:
    workflow.run_pipeline(
        config,
        config_path=config_path,
        jobs=jobs,
        use_snakemake=False if force_python else None,
    )


def _run_validate_data(config: dict) -> None:
    """Validate all datasets in the configuration."""
    failures = []
    successes = []

    raw_dir = Path(config.get("paths", {}).get("raw_dir", "data/raw"))
    ensure_dir(raw_dir)

    logging.info(f"Validating {len(config['datasets'])} datasets...")

    for dataset_info in config["datasets"]:
        dataset_id = dataset_info["id"]
        try:
            logging.info(f"Validating dataset: {dataset_id}")

            # Ensure remote assets are materialized locally if needed
            resolved_entry = download_if_needed(dataset_info, raw_dir, config=config)

            # Check if file exists
            url_path = Path(resolved_entry["url"])
            if not url_path.exists():
                error_msg = f"Dataset file not found: {url_path}"
                logging.error(f"  ✗ {dataset_id}: {error_msg}")
                failures.append((dataset_id, error_msg))
                continue

            # Load and validate the dataset
            adata = load_matrix(resolved_entry, config=config)

            # Validation already happens in load_matrix, but we can do additional checks
            is_valid, error = validate_anndata(
                adata, dataset_id=dataset_id, raise_on_error=False
            )

            if is_valid:
                logging.info(
                    f"  ✓ {dataset_id}: {adata.n_obs} cells, {adata.n_vars} genes"
                )
                successes.append(dataset_id)
            else:
                logging.error(f"  ✗ {dataset_id}: {error}")
                failures.append((dataset_id, error))

        except Exception as e:
            error_msg = str(e)
            logging.error(f"  ✗ {dataset_id}: {error_msg}")
            failures.append((dataset_id, error_msg))

    # Summary
    print("\n" + "=" * 80)
    print(f"Validation Summary: {len(successes)}/{len(config['datasets'])} passed")
    print("=" * 80)

    if successes:
        print(f"\n✓ Passed ({len(successes)}):")
        for dataset_id in successes:
            print(f"  - {dataset_id}")

    if failures:
        print(f"\n✗ Failed ({len(failures)}):")
        for dataset_id, error in failures:
            print(f"  - {dataset_id}")
            print(f"    Error: {error}")

    # Exit with error code if any failures
    if failures:
        sys.exit(1)


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
        "pipeline",
        "validate-data",
        "receptor",
        "tcr",
    ])
    parser.add_argument("--config", default="config/atlas.yaml", help="Config file path")
    parser.add_argument("--dataset", help="Dataset ID (for qc/doublets)")
    parser.add_argument("--log-level", default="INFO", help="Logging level")
    parser.add_argument(
        "--jobs",
        type=int,
        help="Number of parallel jobs when running the Snakemake pipeline",
    )
    parser.add_argument(
        "--force-python",
        action="store_true",
        help="Skip Snakemake and run the Python fallback pipeline",
    )
    parser.add_argument(
        "--stage",
        choices=["ingest", "qc", "analytics", "viz", "all"],
        help="Receptor module stage (receptor command only)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force receptor run even if disabled in config",
    )

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
    elif args.command == "pipeline":
        _run_pipeline_cli(config, args.config, args.jobs, args.force_python)
    elif args.command == "validate-data":
        _run_validate_data(config)
    elif args.command == "receptor":
        receptor_args = argparse.Namespace(
            config=args.config,
            dataset=args.dataset,
            stage=args.stage or "all",
            log_level=args.log_level,
            force=args.force,
        )
        receptor._cli_entry(receptor_args)
    elif args.command == "tcr":
        tcr.run_tcr_analysis(config)
    else:  # pragma: no cover - safeguard
        raise ValueError(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()
