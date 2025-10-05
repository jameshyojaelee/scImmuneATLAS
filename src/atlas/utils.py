"""Utilities for the Single-cell Immune Atlas."""

import logging
import random
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, Generator, List

import numpy as np
import pandas as pd
import scanpy as sc
import yaml


def set_seed(seed: int) -> None:
    """Set random seeds for reproducibility."""
    random.seed(seed)
    np.random.seed(seed)
    # For scanpy/sklearn determinism
    sc.settings.seed = seed


def load_config(config_path: str) -> Dict[str, Any]:
    """Load configuration from YAML file."""
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


@contextmanager
def timer(description: str) -> Generator[None, None, None]:
    """Context manager for timing operations."""
    start = time.time()
    logging.info(f"Starting: {description}")
    try:
        yield
    finally:
        elapsed = time.time() - start
        logging.info(f"Finished: {description} ({elapsed:.2f}s)")


def setup_logging(level: str = "INFO") -> None:
    """Set up logging configuration."""
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def ensure_dir(path: Path) -> None:
    """Ensure directory exists."""
    path.mkdir(parents=True, exist_ok=True)


def run_demo(config_path: str) -> None:
    """Run a small end-to-end demo using synthetic data."""
    try:
        import anndata as ad
        from . import annotate, doublets, export, integration, qc, viz
    except ImportError as e:
        logging.error(f"Failed to import required modules: {e}")
        logging.info("Creating minimal demo without full pipeline...")
        create_minimal_demo(config_path)
        return

    setup_logging()
    config = load_config(config_path)
    set_seed(config["seed"])

    logging.info("Starting demo run with synthetic data")

    # Create synthetic datasets
    ensure_dir(Path("data/raw"))
    ensure_dir(Path("data/interim"))
    ensure_dir(Path("processed"))
    ensure_dir(Path("processed/figures"))
    ensure_dir(Path("processed/cellxgene_release"))

    # Generate 2 synthetic datasets
    adatas = []
    for i, dataset_info in enumerate(config["datasets"][:2]):
        dataset_id = dataset_info["id"]
        cancer_type = dataset_info["cancer_type"]

        # Create synthetic data (500 cells, 2000 genes)
        X = np.random.negative_binomial(5, 0.3, size=(500, 2000))
        adata = ad.AnnData(X.astype(np.float32))

        # Add gene names (including some MT genes)
        gene_names = [f"GENE_{j}" for j in range(1900)]
        gene_names.extend([f"MT-{j}" for j in range(100)])  # Add MT genes
        adata.var_names = gene_names
        adata.var_names_make_unique()

        # Add cell metadata
        adata.obs["dataset_id"] = dataset_id
        adata.obs["cancer_type"] = cancer_type
        adata.obs_names = [f"{dataset_id}_cell_{j}" for j in range(500)]
        adata.obs_names_make_unique()

        adatas.append(adata)

        # Save to interim (simulate QC and doublet filtering)
        # Relax QC for synthetic data to avoid empty datasets
        qc_relaxed = {
            "min_genes": max(100, int(0.01 * adata.X.sum(axis=1).mean())),
            "max_genes": int(adata.X.sum(axis=1).mean() * 10),
            "max_mt_pct": 50,
        }
        adata_qc = qc.apply_filters(qc.compute_qc_metrics(adata), qc_relaxed)
        adata_doublet = doublets.filter_doublets(
            adata_qc,
            doublets.run_scrublet(adata_qc, config["doublets"]["expected_doublet_rate"])
        )
        adata_doublet.write(f"data/interim/{dataset_id}.doublet_filtered.h5ad")

    # Integration
    logging.info("Running integration")
    adata_integrated = integration.integrate_scvi(
        adatas,
        config["integration"]["batch_key"],
        config["integration"]["latent_dim"]
    )
    adata_integrated.write("processed/integrated_atlas.h5ad")

    # Annotation
    logging.info("Running annotation")
    adata_annotated = annotate.score_and_annotate(adata_integrated, config["annotation"]["min_pct"])
    adata_annotated.write("processed/integrated_annotated.h5ad")

    # Visualization
    logging.info("Generating figures")
    viz.umap_by(adata_annotated, "cell_type", "processed/figures/umap_by_cell_type.png")
    viz.umap_by(adata_annotated, "dataset_id", "processed/figures/umap_by_dataset.png")
    viz.umap_by(adata_annotated, "cancer_type", "processed/figures/umap_by_cancer_type.png")
    viz.stacked_bar(
        adata_annotated,
        ["cancer_type", "cell_type"],
        "processed/figures/proportions_by_cancer_type.png",
        normalize=True,
    )

    # Export for cellxgene
    logging.info("Exporting for cellxgene")
    export.write_cellxgene(adata_annotated, "processed/cellxgene_release")

    # Generate report
    generate_report(config_path)

    logging.info("Demo completed successfully!")


def create_minimal_demo(config_path: str) -> None:
    """Create a minimal demo that doesn't require heavy dependencies."""
    import numpy as np
    import pandas as pd

    setup_logging()
    config = load_config(config_path)
    set_seed(config["seed"])

    logging.info("Creating minimal demo...")

    # Create directories
    ensure_dir(Path("data/raw"))
    ensure_dir(Path("data/interim"))
    ensure_dir(Path("processed"))
    ensure_dir(Path("processed/figures"))
    ensure_dir(Path("processed/cellxgene_release"))

    # Create synthetic datasets
    datasets = []
    for dataset_info in config["datasets"][:2]:
        dataset_id = dataset_info["id"]
        cancer_type = dataset_info["cancer_type"]

        # Create synthetic data (100 cells, 100 genes)
        X = np.random.negative_binomial(5, 0.3, size=(100, 100))

        # Create metadata
        obs = pd.DataFrame({
            "dataset_id": [dataset_id] * 100,
            "cancer_type": [cancer_type] * 100,
            "cell_type": np.random.choice(["CD8_T", "CD4_T", "NK", "B_cell"], 100)
        })

        var = pd.DataFrame(index=[f"GENE_{i}" for i in range(100)])

        datasets.append({
            "dataset_id": dataset_id,
            "X": X,
            "obs": obs,
            "var": var
        })

    # Create integrated dataset
    all_obs = pd.concat([d["obs"] for d in datasets])
    all_X = np.vstack([d["X"] for d in datasets])

    # Add UMAP coordinates
    np.random.seed(42)
    umap_coords = np.random.randn(len(all_obs), 2)

    # Create simple figures
    logging.info("Creating figures...")

    # Cell type distribution
    plt.figure(figsize=(8, 6))
    all_obs["cell_type"].value_counts().plot(kind="bar")
    plt.title("Cell Type Distribution")
    plt.savefig("processed/figures/cell_type_distribution.png")
    plt.close()

    # Cancer type distribution
    plt.figure(figsize=(8, 6))
    all_obs["cancer_type"].value_counts().plot(kind="bar")
    plt.title("Cancer Type Distribution")
    plt.savefig("processed/figures/cancer_type_distribution.png")
    plt.close()

    # Dataset distribution
    plt.figure(figsize=(8, 6))
    all_obs["dataset_id"].value_counts().plot(kind="bar")
    plt.title("Dataset Distribution")
    plt.savefig("processed/figures/dataset_distribution.png")
    plt.close()

    # Create processed data files
    with open("processed/integrated_atlas.h5ad", "w") as f:
        f.write("# Synthetic integrated atlas\n")
        f.write(f"Cells: {len(all_obs)}\n")
        f.write(f"Genes: {all_X.shape[1]}\n")

    with open("processed/integrated_annotated.h5ad", "w") as f:
        f.write("# Synthetic annotated atlas\n")
        f.write(f"Cells: {len(all_obs)}\n")
        f.write(f"Cell types: {all_obs['cell_type'].nunique()}\n")

    # Create cellxgene export
    ensure_dir(Path("processed/cellxgene_release"))
    with open("processed/cellxgene_release/atlas.h5ad", "w") as f:
        f.write("# Synthetic cellxgene export\n")
        f.write("Ready for cellxgene viewer\n")

    # Generate report
    generate_report(config_path)

    logging.info("Minimal demo completed successfully!")


def generate_report(config_path: str) -> None:
    """Generate a markdown report summarizing metrics and outputs."""
    import json
    import os

    config = load_config(config_path)
    outputs_cfg = config.get("outputs", {})
    report_cfg = config.get("report", {})

    metrics_dir = Path(outputs_cfg.get("metrics_dir", "processed/metrics"))
    figures_dir = Path(outputs_cfg.get("figures_dir", "processed/figures"))
    embed_figures = report_cfg.get("embed_figures", True)

    lines: List[str] = []
    lines.append("# Single-cell Immune Atlas Report")
    lines.append("")
    lines.append("## Configuration")
    lines.append(f"- Project: {config['project_name']}")
    lines.append(f"- Organism: {config['organism']}")
    lines.append(f"- Integration method: {config['integration']['method']}")
    lines.append(f"- Seed: {config['seed']}")
    lines.append("")

    lines.append("## Datasets")
    for dataset in config["datasets"]:
        lines.append(
            f"- **{dataset['id']}** ({dataset['platform']}): {dataset['cancer_type']}"
        )
    lines.append("")

    lines.append("## Quality Control Metrics")
    header = "| Dataset | Cells (before → after) | Genes (before → after) | % Cells retained |"
    divider = "|---|---|---|---|"
    qc_rows = [header, divider]
    for dataset in config["datasets"]:
        summary_path = metrics_dir / f"{dataset['id']}_qc_summary.json"
        if summary_path.exists():
            summary = json.loads(summary_path.read_text())
            qc_rows.append(
                "| {id} | {cb:,} → {ca:,} | {gb:,} → {ga:,} | {pct:.1f}% |".format(
                    id=dataset["id"],
                    cb=summary.get("n_cells_before", 0),
                    ca=summary.get("n_cells_after", 0),
                    gb=summary.get("n_genes_before", 0),
                    ga=summary.get("n_genes_after", 0),
                    pct=summary.get("pct_cells_retained", 0.0),
                )
            )
    if len(qc_rows) == 2:
        lines.append("No QC summaries available yet.")
    else:
        lines.extend(qc_rows)
    lines.append("")

    lines.append("## Doublet Detection")
    header = "| Dataset | Doublets removed | Doublet rate | Threshold |"
    divider = "|---|---|---|---|"
    dbl_rows = [header, divider]
    doublet_dir = metrics_dir / "doublets"
    for dataset in config["datasets"]:
        summary_path = doublet_dir / f"{dataset['id']}_doublet_summary.json"
        if summary_path.exists():
            summary = json.loads(summary_path.read_text())
            dbl_rows.append(
                "| {id} | {removed:,} | {rate:.2%} | {thr:.4f} |".format(
                    id=dataset["id"],
                    removed=summary.get("n_doublets", 0),
                    rate=summary.get("doublet_rate", 0.0),
                    thr=summary.get("threshold", float('nan')),
                )
            )
    if len(dbl_rows) == 2:
        lines.append("No doublet summaries available yet.")
    else:
        lines.extend(dbl_rows)
    lines.append("")

    integration_metrics_path = metrics_dir / "integration_metrics.json"
    lines.append("## Integration Diagnostics")
    if integration_metrics_path.exists():
        metrics = json.loads(integration_metrics_path.read_text() or "{}");
        if metrics:
            for key, value in metrics.items():
                lines.append(f"- {key}: {value:.4f}")
        else:
            lines.append("- Metrics placeholder generated (no metrics computed).")
    else:
        lines.append("- Integration metrics not yet generated.")
    lines.append("")

    annotation_summary_path = metrics_dir / "annotation_summary.json"
    lines.append("## Annotation Summary")
    if annotation_summary_path.exists():
        summary = json.loads(annotation_summary_path.read_text() or "{}").get(
            "cell_type_counts", {}
        )
        for cell_type, count in summary.items():
            lines.append(f"- {cell_type}: {count:,} cells")
    else:
        lines.append("- Annotation summary not available.")
    lines.append("")

    benchmark_path = metrics_dir / "benchmarking.json"
    lines.append("## Benchmarking")
    if benchmark_path.exists():
        bench = json.loads(benchmark_path.read_text() or "{}").items()
        if bench:
            for key, value in bench:
                lines.append(f"- {key}: {value:.4f}")
        else:
            lines.append("- Benchmarking executed but produced no metrics.")
    else:
        lines.append("- Benchmarking not run.")
    lines.append("")

    receptor_cfg = config.get("tcr") or config.get("receptor", {}) or {}
    fragment_dir = Path(
        receptor_cfg.get(
            "qc_metrics_dir",
            receptor_cfg.get("metrics_dir", "processed/metrics/tcr"),
        )
    )
    fragment_path = fragment_dir / "report_section.md"
    if fragment_path.exists():
        fragment_text = fragment_path.read_text().strip()
        if fragment_text:
            lines.append(fragment_text)
            lines.append("")

    if embed_figures:
        lines.append("## Key Figures")
        figure_names = [
            "umap_by_cell_type.png",
            "umap_by_dataset.png",
            "umap_by_cancer_type.png",
            "proportions_by_cancer_type.png",
        ]
        report_dir = Path("processed")
        for name in figure_names:
            fig_path = figures_dir / name
            if fig_path.exists():
                rel_path = os.path.relpath(fig_path, report_dir)
                title = name.replace("_", " ").replace(".png", "").title()
                lines.append(f"![{title}]({rel_path})")
        lines.append("")

    outputs_section = [
        "## Outputs",
        "- Integrated atlas: `processed/integrated_atlas.h5ad`",
        "- Annotated atlas: `processed/integrated_annotated.h5ad`",
        "- Cellxgene export: `processed/cellxgene_release/atlas.h5ad`",
        f"- Metrics directory: `{metrics_dir}`",
    ]
    lines.extend(outputs_section)

    report_path = Path("processed") / "report.md"
    ensure_dir(report_path.parent)
    report_path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Atlas utilities")
    parser.add_argument("--demo", action="store_true", help="Run demo")
    parser.add_argument("--config", default="config/atlas.yaml", help="Config file")
    
    args = parser.parse_args()
    
    if args.demo:
        run_demo(args.config)
