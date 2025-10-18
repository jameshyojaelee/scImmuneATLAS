"""Smoke tests for demo utilities."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

anndata = pytest.importorskip("anndata")  # noqa: F401  (ensures availability for demo tests)
import anndata as ad  # type: ignore  # noqa: E402

from src.atlas.utils import create_minimal_demo, run_demo


def _demo_config_dict() -> dict:
    """Return a lightweight configuration for demo tests."""
    return {
        "project_name": "demo",
        "seed": 42,
        "organism": "human",
        "datasets": [
            {
                "id": "DEMO_A",
                "modality": "scRNA",
                "url": "",
                "platform": "synthetic",
                "cancer_type": "Melanoma",
            },
            {
                "id": "DEMO_B",
                "modality": "scRNA",
                "url": "",
                "platform": "synthetic",
                "cancer_type": "NSCLC",
            },
        ],
        "qc": {
            "min_genes": 200,
            "max_genes": 6000,
            "max_mt_pct": 20,
            "min_counts": 500,
            "max_counts": 50000,
            "min_cells_per_gene": 3,
            "qc_plots_dir": "processed/figures/qc",
        },
        "doublets": {
            "method": "scrublet",
            "expected_doublet_rate": 0.06,
            "mask_high_mt_genes": True,
            "save_plots": True,
            "scrublet_kwargs": {
                "min_counts": 2,
                "min_cells": 3,
                "min_gene_variability_pctl": 85,
                "n_prin_comps": 30,
            },
        },
        "integration": {
            "method": "harmony",
            "latent_dim": 20,
            "batch_key": "dataset_id",
            "max_epochs": 25,
            "use_gpu": "auto",
            "shared_hvg": True,
            "n_hvg": 800,
            "metrics": ["silhouette"],
        },
        "neighbors": {"n_neighbors": 15, "n_pcs": 30},
        "umap": {"min_dist": 0.4, "spread": 1.0},
        "annotation": {
            "strategy": "markers",
            "min_pct": 0.05,
            "score_method": "scanpy_score_genes",
            "min_confidence_pct": 0.7,
            "marker_genes": "src/atlas/markers/immune_markers_human.tsv",
        },
        "outputs": {
            "cellxgene_dir": "processed/cellxgene_release",
            "figures_dir": "processed/figures",
            "metrics_dir": "processed/metrics",
        },
        "benchmarking": {"enabled": False},
        "tcr": {"enabled": False},
    }


@pytest.fixture()
def demo_config_path(tmp_path: Path) -> Path:
    config_path = tmp_path / "demo_config.yaml"
    config_path.write_text(yaml.safe_dump(_demo_config_dict()))
    return config_path


def _assert_common_artifacts(base_path: Path) -> None:
    assert (base_path / "processed/integrated_atlas.h5ad").exists()
    assert (base_path / "processed/integrated_annotated.h5ad").exists()
    assert (base_path / "processed/cellxgene_release/atlas.h5ad").exists()
    assert (base_path / "processed/report.md").exists()


def test_run_demo_executes_pipeline(monkeypatch: pytest.MonkeyPatch, demo_config_path: Path, tmp_path: Path) -> None:
    monkeypatch.chdir(tmp_path)
    run_demo(str(demo_config_path))

    _assert_common_artifacts(tmp_path)

    metrics_path = tmp_path / "processed/metrics/integration_metrics.json"
    assert metrics_path.exists()
    metrics = json.loads(metrics_path.read_text())
    assert "method" in metrics

    integrated = ad.read_h5ad(tmp_path / "processed/integrated_atlas.h5ad")
    assert integrated.n_obs > 0
    assert integrated.obs["dataset_id"].nunique() == 2


def test_create_minimal_demo_generates_h5ad(monkeypatch: pytest.MonkeyPatch, demo_config_path: Path, tmp_path: Path) -> None:
    monkeypatch.chdir(tmp_path)
    create_minimal_demo(str(demo_config_path))

    _assert_common_artifacts(tmp_path)

    atlas = ad.read_h5ad(tmp_path / "processed/integrated_atlas.h5ad")
    assert atlas.n_obs == 80
