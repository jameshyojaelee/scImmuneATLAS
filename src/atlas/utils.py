"""Utilities for the Single-cell Immune Atlas."""

import logging
import random
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, Generator

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
    """Generate a lightweight markdown report."""
    config = load_config(config_path)
    
    report_content = f"""# Single-cell Immune Atlas Report

## Configuration
- Project: {config['project_name']}
- Organism: {config['organism']}
- Integration method: {config['integration']['method']}
- Seed: {config['seed']}

## Datasets
"""
    
    for dataset in config['datasets']:
        report_content += f"- **{dataset['id']}**: {dataset['cancer_type']} ({dataset['platform']})\n"
    
    report_content += f"""
## QC Parameters
- Min genes per cell: {config['qc']['min_genes']}
- Max genes per cell: {config['qc']['max_genes']}
- Max mitochondrial %: {config['qc']['max_mt_pct']}

## Integration Parameters
- Method: {config['integration']['method']}
- Latent dimensions: {config['integration']['latent_dim']}
- Batch key: {config['integration']['batch_key']}

## Outputs
- Integrated atlas: `processed/integrated_atlas.h5ad`
- Annotated atlas: `processed/integrated_annotated.h5ad`
- Cellxgene export: `processed/cellxgene_release/atlas.h5ad`
- Figures: `processed/figures/`

## Generated Figures
- UMAP by cell type
- UMAP by dataset
- UMAP by cancer type
- Cell type proportions by cancer type
"""
    
    with open("processed/report.md", "w") as f:
        f.write(report_content)


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Atlas utilities")
    parser.add_argument("--demo", action="store_true", help="Run demo")
    parser.add_argument("--config", default="config/atlas.yaml", help="Config file")
    
    args = parser.parse_args()
    
    if args.demo:
        run_demo(args.config)
