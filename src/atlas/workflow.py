"""Workflow orchestration helpers for scImmuneATLAS."""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

from . import annotate, benchmark, doublets, export, integration, qc, viz, utils, receptor, tcr
from .receptor.config import global_receptor_config
from .utils import ensure_dir, timer


def dataset_ids(config: Dict) -> List[str]:
    """Return dataset identifiers declared in the configuration."""
    return [dataset["id"] for dataset in config.get("datasets", [])]


def collect_pipeline_targets(config: Dict) -> List[str]:
    """Materialized artifacts that define a successful pipeline run."""
    outputs_cfg = config.get("outputs", {})
    metrics_dir = Path(outputs_cfg.get("metrics_dir", "processed/metrics"))
    figures_dir = Path(outputs_cfg.get("figures_dir", "processed/figures"))
    cellxgene_dir = Path(outputs_cfg.get("cellxgene_dir", "processed/cellxgene_release"))

    targets: List[Path] = [
        Path("processed/integrated_annotated.h5ad"),
        cellxgene_dir / "atlas.h5ad",
        figures_dir / "umap_by_cell_type.png",
        figures_dir / "umap_by_dataset.png",
        figures_dir / "umap_by_cancer_type.png",
        figures_dir / "proportions_by_cancer_type.png",
        Path("processed/report.md"),
        metrics_dir / "integration_metrics.json",
        metrics_dir / "annotation_summary.json",
    ]

    receptor_cfg = global_receptor_config(config)
    if receptor_cfg.get("enabled"):
        receptor_metrics_dir = Path(
            receptor_cfg.get("qc_metrics_dir", receptor_cfg.get("metrics_dir", "processed/metrics/tcr"))
        )
        receptor_figures_dir = Path(
            receptor_cfg.get("figures_dir", "processed/figures/tcr")
        )
        targets.extend(
            [
                receptor_metrics_dir / "repertoire_summary.json",
                receptor_metrics_dir / "clonotype_expansion.tsv",
                receptor_metrics_dir / "diversity_indices.json",
                receptor_metrics_dir / "vj_usage.tsv",
                receptor_metrics_dir / "chain_pairing_summary.json",
                receptor_metrics_dir / "v_usage_summary.json",
                receptor_metrics_dir / "report_section.md",
                receptor_figures_dir / "clonotype_frequency.png",
                receptor_figures_dir / "umap_clonal_expansion.png",
                receptor_figures_dir / "cdr3_spectratype.png",
                receptor_figures_dir / "vj_usage_heatmap.png",
            ]
        )

    if receptor_cfg.get("enabled"):
        tcr_metrics_dir = Path(receptor_cfg.get("qc_metrics_dir", "processed/metrics/tcr"))
        targets.append(tcr_metrics_dir / "tcr_summary.json")

    if config.get("benchmarking", {}).get("enabled", False):
        targets.append(metrics_dir / "benchmarking.json")

    return [str(path) for path in targets]


def _attach_file_logger(log_dir: Path, filename: str = "pipeline.log") -> logging.Handler:
    """Attach a file handler to the root logger and return it."""
    ensure_dir(log_dir)
    handler = logging.FileHandler(log_dir / filename)
    handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logging.getLogger().addHandler(handler)
    return handler


def _remove_handler(handler: logging.Handler) -> None:
    handler.flush()
    handler.close()
    logging.getLogger().removeHandler(handler)


def run_pipeline(
    config: Dict,
    *,
    config_path: str = "config/atlas.yaml",
    use_snakemake: Optional[bool] = None,
    jobs: Optional[int] = None,
    log_dir: Optional[Path | str] = "logs",
) -> None:
    """Execute the full atlas workflow.

    If Snakemake is available (or explicitly requested) the workflow is delegated to
    the Snakefile. Otherwise, the individual Python stages are executed sequentially.
    """

    logging.info("Starting pipeline orchestration")

    log_path = Path(log_dir) if log_dir else Path("logs")
    file_handler = _attach_file_logger(log_path)
    try:
        if use_snakemake is None:
            use_snakemake = shutil.which("snakemake") is not None

        if use_snakemake:
            cmd = ["snakemake"]
            if jobs:
                cmd += ["-j", str(jobs)]
            logging.info("Running Snakemake command: %s", " ".join(cmd))
            subprocess.run(cmd, check=True)
            logging.info("Snakemake workflow finished successfully")
            return

        logging.info("Snakemake not available; executing Python fallback pipeline")
        receptor_enabled = receptor.is_enabled(config)
        receptor_datasets = set(receptor.datasets_with_receptor(config)) if receptor_enabled else set()
        for dataset_id in dataset_ids(config):
            with timer(f"QC ({dataset_id})"):
                qc.process_dataset_qc(dataset_id, config)
            with timer(f"Doublets ({dataset_id})"):
                doublets.process_dataset_doublets(dataset_id, config)
            if receptor_enabled and dataset_id in receptor_datasets:
                with timer(f"Receptor ingest ({dataset_id})"):
                    receptor.ingest_dataset(dataset_id, config)
                with timer(f"Receptor QC ({dataset_id})"):
                    receptor.qc_dataset(dataset_id, config)

        with timer("Integration"):
            integration.run_integration(config)
        with timer("Annotation"):
            annotate.run_annotation(config)
        if receptor_enabled:
            with timer("Receptor analytics"):
                receptor.run_all(config)
            with timer("TCR analysis"):
                tcr.run_tcr_analysis(config)
        with timer("Export"):
            export.run_export(config)
        with timer("Visualization"):
            viz.generate_all_figures(config)

        if config.get("benchmarking", {}).get("enabled", False):
            with timer("Benchmarking"):
                benchmark.run_benchmark(config)

        with timer("Report"):
            utils.generate_report(config_path)

        logging.info("Python fallback pipeline completed successfully")

    finally:
        _remove_handler(file_handler)
