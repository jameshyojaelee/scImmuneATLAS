"""Benchmarking utilities for the Single-cell Immune Atlas."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict

import anndata as ad
import numpy as np

from .utils import ensure_dir, load_config, setup_logging, timer


def _compute_label_metrics(
    adata: ad.AnnData,
    ref_adata: ad.AnnData,
    label_key: str,
    metrics: Dict[str, float],
) -> None:
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    common = np.intersect1d(adata.obs_names, ref_adata.obs_names)
    if common.size == 0:
        logging.warning("No overlapping cells between atlas and reference; skipping label metrics")
        return

    atlas_labels = adata.obs.loc[common, "cell_type"].astype(str)
    reference_labels = ref_adata.obs.loc[common, label_key].astype(str)

    metrics["ari"] = float(adjusted_rand_score(reference_labels, atlas_labels))
    metrics["nmi"] = float(normalized_mutual_info_score(reference_labels, atlas_labels))


def run_benchmark(config: Dict) -> None:
    benchmarking_cfg = config.get("benchmarking", {})
    if not benchmarking_cfg.get("enabled", False):
        logging.info("Benchmarking disabled; skipping benchmark step")
        return

    integrated_path = Path("processed") / "integrated_annotated.h5ad"
    if not integrated_path.exists():
        raise FileNotFoundError(
            "Integrated annotated atlas not found. Run annotation before benchmarking."
        )

    adata = ad.read_h5ad(integrated_path)
    metrics: Dict[str, float] = {}

    reference_path = benchmarking_cfg.get("reference_atlas")
    if reference_path:
        ref_path = Path(reference_path)
        if ref_path.exists():
            ref_adata = ad.read_h5ad(ref_path)
            reference_label = benchmarking_cfg.get("reference_label", "cell_type")
            if reference_label not in ref_adata.obs.columns:
                logging.warning(
                    "Reference label column '%s' missing; available columns: %s",
                    reference_label,
                    list(ref_adata.obs.columns),
                )
            else:
                _compute_label_metrics(adata, ref_adata, reference_label, metrics)
        else:
            logging.warning("Reference atlas %s not found; skipping label benchmarking", ref_path)

    metrics_dir = Path(config["outputs"].get("metrics_dir", "processed/metrics"))
    ensure_dir(metrics_dir)
    benchmark_path = metrics_dir / "benchmarking.json"
    with open(benchmark_path, "w") as f:
        json.dump(metrics, f, indent=2)
    logging.info("Saved benchmarking results to %s", benchmark_path)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run benchmarking against reference atlases")
    parser.add_argument("--config", default="config/atlas.yaml", help="Config file path")

    args = parser.parse_args()

    setup_logging()
    cfg = load_config(args.config)

    with timer("Benchmarking"):
        run_benchmark(cfg)
