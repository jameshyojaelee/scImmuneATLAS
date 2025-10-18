"""Lightweight regression tests for pipeline orchestration."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Dict

import numpy as np
import pytest
import yaml
import sys

import anndata as ad  # type: ignore

import src.atlas.workflow as workflow


@pytest.fixture()
def synthetic_pipeline_assets(tmp_path: Path) -> Dict[str, Path]:
    """Create synthetic datasets and a minimal config pointing at them."""
    raw_dir = tmp_path / "data" / "raw"
    interim_dir = tmp_path / "data" / "interim"
    processed_dir = tmp_path / "processed"
    logs_dir = tmp_path / "logs"
    metrics_dir = processed_dir / "metrics"
    figures_dir = processed_dir / "figures"
    cellxgene_dir = processed_dir / "cellxgene_release"

    for directory in (raw_dir, interim_dir, processed_dir, metrics_dir, figures_dir, cellxgene_dir, logs_dir):
        directory.mkdir(parents=True, exist_ok=True)

    def _write_dataset(dataset_id: str, seed: int) -> Path:
        rng = np.random.default_rng(seed)
        X = rng.poisson(lam=2.0, size=(25, 40)).astype(np.float32)
        adata = ad.AnnData(X)
        adata.var_names = [f"GENE_{i}" for i in range(40)]
        adata.obs_names = [f"{dataset_id}_cell_{i}" for i in range(25)]
        adata.obs["dataset_id"] = dataset_id
        adata.obs["cancer_type"] = f"Cancer_{seed}"
        adata.obs["platform"] = "synthetic"
        path = raw_dir / f"{dataset_id.lower()}.h5ad"
        adata.write_h5ad(path)
        return path

    dataset_a = _write_dataset("DATASET_A", seed=1)
    dataset_b = _write_dataset("DATASET_B", seed=2)

    config = {
        "project_name": "workflow-regression",
        "seed": 123,
        "organism": "human",
        "datasets": [
            {
                "id": "DATASET_A",
                "modality": "scRNA",
                "url": str(dataset_a),
                "platform": "synthetic",
                "cancer_type": "Melanoma",
            },
            {
                "id": "DATASET_B",
                "modality": "scRNA",
                "url": str(dataset_b),
                "platform": "synthetic",
                "cancer_type": "NSCLC",
            },
        ],
        "qc": {
            "min_genes": 50,
            "max_genes": 5000,
            "max_mt_pct": 40,
            "min_counts": 0,
            "max_counts": 100000,
            "min_cells_per_gene": 1,
            "qc_plots_dir": str(figures_dir / "qc"),
        },
        "doublets": {
            "method": "scrublet",
            "expected_doublet_rate": 0.04,
            "mask_high_mt_genes": True,
            "save_plots": True,
            "scrublet_kwargs": {
                "min_counts": 2,
                "min_cells": 3,
                "min_gene_variability_pctl": 85,
                "n_prin_comps": 20,
            },
        },
        "integration": {
            "method": "harmony",
            "latent_dim": 15,
            "batch_key": "dataset_id",
            "shared_hvg": True,
            "n_hvg": 800,
            "metrics": ["silhouette"],
        },
        "neighbors": {"n_neighbors": 10, "n_pcs": 20},
        "umap": {"min_dist": 0.5, "spread": 1.0},
        "annotation": {
            "strategy": "markers",
            "min_pct": 0.05,
            "score_method": "scanpy_score_genes",
            "min_confidence_pct": 0.6,
            "marker_genes": "src/atlas/markers/immune_markers_human.tsv",
        },
        "outputs": {
            "metrics_dir": str(metrics_dir),
            "figures_dir": str(figures_dir),
            "cellxgene_dir": str(cellxgene_dir),
        },
        "benchmarking": {"enabled": False},
        "tcr": {"enabled": False},
    }

    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config))

    return {
        "config_path": config_path,
        "base_dir": tmp_path,
        "metrics_dir": metrics_dir,
        "figures_dir": figures_dir,
        "processed_dir": processed_dir,
        "logs_dir": logs_dir,
    }


@pytest.fixture()
def patched_workflow(monkeypatch: pytest.MonkeyPatch):
    """Patch expensive workflow steps to keep tests lightweight."""

    def _noop(*args, **kwargs):
        return None

    def _write_stub_metrics(config: Dict):
        metrics_dir = Path(config["outputs"]["metrics_dir"])
        metrics_dir.mkdir(parents=True, exist_ok=True)
        (metrics_dir / "integration_metrics.json").write_text("method: stub\n")

    def _write_stub_annotation(config: Dict):
        processed = Path("processed")
        processed.mkdir(parents=True, exist_ok=True)
        atlas_path = processed / "integrated_atlas.h5ad"
        if not atlas_path.exists():
            adata = ad.AnnData(np.eye(5, dtype=np.float32))
            adata.write_h5ad(atlas_path)
        annotated_path = processed / "integrated_annotated.h5ad"
        ad.read_h5ad(atlas_path).write_h5ad(annotated_path)
        metrics_dir = Path(config["outputs"]["metrics_dir"])
        metrics_dir.mkdir(parents=True, exist_ok=True)
        (metrics_dir / "annotation_summary.json").write_text("cell_type_counts: {}\n")

    monkeypatch.setattr(workflow.integration, "run_integration", _write_stub_metrics)
    monkeypatch.setattr(workflow.annotate, "run_annotation", _write_stub_annotation)
    monkeypatch.setattr(workflow.export, "run_export", _noop)
    monkeypatch.setattr(workflow.viz, "generate_all_figures", _noop)


def _assert_pipeline_outputs(base_dir: Path) -> None:
    processed = base_dir / "processed"
    metrics = processed / "metrics"
    figures = processed / "figures"
    logs = base_dir / "logs" / "pipeline.log"

    assert (processed / "integrated_atlas.h5ad").exists()
    assert (processed / "integrated_annotated.h5ad").exists()
    assert (metrics / "integration_metrics.json").exists()
    assert (metrics / "annotation_summary.json").exists()
    assert logs.exists()


def test_python_fallback_pipeline(monkeypatch: pytest.MonkeyPatch, synthetic_pipeline_assets, patched_workflow):
    monkeypatch.chdir(synthetic_pipeline_assets["base_dir"])
    config_path = synthetic_pipeline_assets["config_path"]
    cfg = yaml.safe_load(config_path.read_text())

    workflow.run_pipeline(
        cfg,
        config_path=str(synthetic_pipeline_assets["config_path"]),
        use_snakemake=False,
        jobs=1,
        log_dir=synthetic_pipeline_assets["logs_dir"],
    )
    _assert_pipeline_outputs(synthetic_pipeline_assets["base_dir"])

    cfg = yaml.safe_load(config_path.read_text())
    workflow.run_pipeline(
        cfg,
        config_path=str(synthetic_pipeline_assets["config_path"]),
        use_snakemake=False,
        jobs=1,
        log_dir=synthetic_pipeline_assets["logs_dir"],
    )
    _assert_pipeline_outputs(synthetic_pipeline_assets["base_dir"])


def test_cli_pipeline_force_python(monkeypatch: pytest.MonkeyPatch, synthetic_pipeline_assets, patched_workflow):
    monkeypatch.chdir(synthetic_pipeline_assets["base_dir"])
    config_path = synthetic_pipeline_assets["config_path"]
    subprocess.check_call(
        [
            sys.executable,
            "-m",
            "src.atlas.cli",
            "pipeline",
            "--config",
            str(config_path),
            "--force-python",
        ]
    )
    _assert_pipeline_outputs(synthetic_pipeline_assets["base_dir"])
