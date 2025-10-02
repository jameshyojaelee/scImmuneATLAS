"""Tests for workflow orchestration."""

from pathlib import Path

import src.atlas.workflow as workflow


def test_run_pipeline_python_fallback(monkeypatch, tmp_path):
    calls = []

    def record(stage):
        def _inner(*args, **kwargs):
            calls.append(stage)
        return _inner

    def record_with_dataset(stage):
        def _inner(dataset_id, config):
            calls.append((stage, dataset_id))
        return _inner

    monkeypatch.setattr(workflow.qc, "process_dataset_qc", record_with_dataset("qc"))
    monkeypatch.setattr(workflow.doublets, "process_dataset_doublets", record_with_dataset("doublets"))
    monkeypatch.setattr(workflow.integration, "run_integration", record("integration"))
    monkeypatch.setattr(workflow.annotate, "run_annotation", record("annotation"))
    monkeypatch.setattr(workflow.export, "run_export", record("export"))
    monkeypatch.setattr(workflow.viz, "generate_all_figures", record("viz"))
    monkeypatch.setattr(workflow.benchmark, "run_benchmark", record("benchmark"))
    monkeypatch.setattr(
        workflow.utils,
        "generate_report",
        lambda config_path: calls.append(("report", config_path)),
    )

    config = {
        "datasets": [{"id": "A"}, {"id": "B"}],
        "benchmarking": {"enabled": True},
    }

    log_dir = tmp_path / "logs"
    workflow.run_pipeline(
        config,
        config_path="config/atlas.yaml",
        use_snakemake=False,
        log_dir=log_dir,
    )

    assert calls == [
        ("qc", "A"),
        ("doublets", "A"),
        ("qc", "B"),
        ("doublets", "B"),
        "integration",
        "annotation",
        "export",
        "viz",
        "benchmark",
        ("report", "config/atlas.yaml"),
    ]

    log_file = Path(log_dir) / "pipeline.log"
    assert log_file.exists()


def test_run_pipeline_snakemake(monkeypatch, tmp_path):
    recorded = {}

    def fake_run(cmd, check):
        recorded["cmd"] = cmd
        recorded["check"] = check

    monkeypatch.setattr(workflow.subprocess, "run", fake_run)

    config = {"datasets": []}

    workflow.run_pipeline(
        config,
        config_path="config/atlas.yaml",
        use_snakemake=True,
        jobs=4,
        log_dir=tmp_path / "logs",
    )

    assert recorded["cmd"] == ["snakemake", "-j", "4"]
    assert recorded["check"] is True
*** End Patch
