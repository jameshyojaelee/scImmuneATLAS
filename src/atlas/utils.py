"""Utilities for the Single-cell Immune Atlas."""

import copy
import json
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


def _ensure_demo_outputs(config: Dict[str, Any]) -> None:
    ensure_dir(Path("data/raw"))
    ensure_dir(Path("data/interim"))
    ensure_dir(Path("processed"))

    outputs = config.get("outputs", {})
    ensure_dir(Path(outputs.get("metrics_dir", "processed/metrics")))
    ensure_dir(Path(outputs.get("figures_dir", "processed/figures")))
    ensure_dir(Path(outputs.get("cellxgene_dir", "processed/cellxgene_release")))


def _write_demo_config_file(config: Dict[str, Any]) -> Path:
    path = Path("processed") / "demo_config.yaml"
    ensure_dir(path.parent)
    path.write_text(yaml.safe_dump(config, sort_keys=False))
    return path


def _synthesize_demo_dataset(
    dataset_id: str,
    cancer_type: str,
    *,
    rng: np.random.Generator,
    n_cells: int,
    n_genes: int,
):
    import anndata as ad

    marker_panel = {
        "CD8_T": ["CD8A", "NKG7"],
        "CD4_T": ["CD4", "IL7R"],
        "Treg": ["FOXP3", "IL2RA"],
        "NK": ["GNLY", "PRF1"],
        "B_cell": ["MS4A1", "CD79A"],
        "Plasma": ["MZB1", "SDC1"],
        "Mono": ["LST1", "LYZ"],
        "DC_cDC1": ["CLEC9A"],
        "DC_cDC2": ["FCER1A", "ITGAX"],
    }
    marker_genes = sorted({gene for genes in marker_panel.values() for gene in genes})
    mt_genes = [f"MT-{i}" for i in range(max(10, min(40, n_genes // 20)))]

    remaining_genes = n_genes - len(marker_genes) - len(mt_genes)
    extra_genes = [f"GENE_{i}" for i in range(max(0, remaining_genes))]
    var_names = (marker_genes + mt_genes + extra_genes)[:n_genes]
    var_name_to_idx = {name: idx for idx, name in enumerate(var_names)}

    counts = rng.poisson(lam=1.5, size=(n_cells, n_genes)).astype(np.float32)
    counts += rng.binomial(1, 0.35, size=(n_cells, n_genes)).astype(np.float32)

    cell_types = rng.choice(list(marker_panel.keys()), size=n_cells, replace=True)
    for idx, cell_type in enumerate(cell_types):
        expressed = rng.choice(n_genes, size=rng.integers(low=350, high=min(n_genes, 600)), replace=False)
        counts[idx, expressed] += rng.integers(1, 6, size=expressed.size)
        for gene in marker_panel[cell_type]:
            gene_idx = var_name_to_idx.get(gene)
            if gene_idx is not None:
                counts[idx, gene_idx] += rng.integers(6, 12)

    adata = ad.AnnData(counts)
    adata.var_names = var_names
    adata.var_names_make_unique()
    adata.var["feature_name"] = adata.var_names

    adata.obs_names = [f"{dataset_id}_cell_{i}" for i in range(n_cells)]
    adata.obs_names_make_unique()
    adata.obs["dataset_id"] = dataset_id
    adata.obs["cancer_type"] = cancer_type
    adata.obs["platform"] = "synthetic"
    adata.obs["cell_type"] = cell_types

    return adata


def _run_fallback_integration(config: Dict[str, Any], dataset_ids: List[str]) -> None:
    import anndata as ad

    logging.info("Executing fallback PCA/UMAP integration for demo data")
    adatas = []
    for dataset_id in dataset_ids:
        input_path = Path("data/interim") / f"{dataset_id}.doublet_filtered.h5ad"
        if not input_path.exists():
            raise FileNotFoundError(f"Expected doublet-filtered dataset missing: {input_path}")
        adatas.append(ad.read_h5ad(input_path))

    if not adatas:
        raise ValueError("Fallback integration received no datasets")

    combined = ad.concat(adatas, join="outer", label="dataset_id", keys=dataset_ids, index_unique="-")
    combined.var_names_make_unique()

    seed = config.get("seed", 0)
    batch_key = config.get("integration", {}).get("batch_key", "dataset_id")
    if batch_key not in combined.obs.columns:
        combined.obs[batch_key] = combined.obs.get("dataset_id")

    with timer("Fallback PCA/UMAP integration"):
        sc.pp.normalize_total(combined, target_sum=1e4)
        sc.pp.log1p(combined)
        try:
            sc.pp.highly_variable_genes(
                combined,
                n_top_genes=min(combined.n_vars, 2000),
                batch_key=batch_key,
            )
            hv_mask = combined.var.get("highly_variable")
            if hv_mask is not None and hv_mask.sum() >= 10:
                combined = combined[:, hv_mask].copy()
        except Exception as exc:  # pragma: no cover - defensive
            logging.warning("Fallback HVG selection failed; proceeding without HVG filtering (%s)", exc)
        sc.pp.scale(combined, max_value=10)
        sc.tl.pca(combined, svd_solver="arpack", random_state=seed)
        sc.pp.neighbors(
            combined,
            n_neighbors=config.get("neighbors", {}).get("n_neighbors", 15),
            random_state=seed,
        )
        sc.tl.umap(
            combined,
            min_dist=config.get("umap", {}).get("min_dist", 0.4),
            spread=config.get("umap", {}).get("spread", 1.0),
            random_state=seed,
        )

    integrated_path = Path("processed") / "integrated_atlas.h5ad"
    ensure_dir(integrated_path.parent)
    combined.write(integrated_path)

    metrics_dir = Path(config.get("outputs", {}).get("metrics_dir", "processed/metrics"))
    ensure_dir(metrics_dir)
    fallback_metrics = {
        "method": "fallback",
        "details": "Fallback PCA/UMAP integration executed (demo)",
    }
    (metrics_dir / "integration_metrics.json").write_text(json.dumps(fallback_metrics, indent=2))
    config.setdefault("integration", {})["method"] = "fallback"


def run_demo(config_path: str, *, max_datasets: int = 2) -> None:
    """Run an end-to-end demo using synthetic data and the standard pipeline."""
    try:
        import anndata  # noqa: F401
        from . import annotate, doublets, export, integration, qc, viz
    except ImportError as exc:
        logging.error("Demo run requires optional dependencies: %s", exc)
        logging.info("Falling back to minimal demo assets.")
        create_minimal_demo(config_path)
        return

    setup_logging()
    original_config = load_config(config_path)
    demo_config = copy.deepcopy(original_config)

    seed = demo_config.get("seed", 0)
    set_seed(seed)

    dataset_entries = list(demo_config.get("datasets", []))[:max_datasets]
    if not dataset_entries:
        raise ValueError("Configuration must contain at least one dataset for the demo.")

    demo_config["datasets"] = dataset_entries
    demo_config.setdefault("outputs", {})
    demo_config["outputs"].setdefault("metrics_dir", "processed/metrics")
    demo_config["outputs"].setdefault("figures_dir", "processed/figures")
    demo_config["outputs"].setdefault("cellxgene_dir", "processed/cellxgene_release")
    demo_config.setdefault("benchmarking", {"enabled": False})
    demo_config.setdefault("tcr", {"enabled": False})

    _ensure_demo_outputs(demo_config)

    rng = np.random.default_rng(seed)
    n_cells = max(160, 80 * len(dataset_entries))
    n_genes = 900

    for dataset in dataset_entries:
        dataset_id = dataset["id"]
        cancer_type = dataset.get("cancer_type", "Unknown")
        adata = _synthesize_demo_dataset(
            dataset_id,
            cancer_type,
            rng=rng,
            n_cells=n_cells,
            n_genes=n_genes,
        )
        out_path = Path("data/raw") / f"{dataset_id.lower()}_demo.h5ad"
        adata.write_h5ad(out_path)
        dataset["url"] = str(out_path)
        dataset["platform"] = dataset.get("platform", "synthetic")
        for key in ("receptor_path", "receptor_format", "receptor_sha256", "receptor"):
            dataset.pop(key, None)

    qc_cfg = copy.deepcopy(demo_config.get("qc", {}))
    qc_cfg["min_genes"] = min(qc_cfg.get("min_genes", 200), max(150, n_genes // 3))
    qc_cfg["max_genes"] = max(qc_cfg.get("max_genes", 6000), qc_cfg["min_genes"] + 100)
    qc_cfg["max_mt_pct"] = max(qc_cfg.get("max_mt_pct", 15), 30)
    qc_cfg.setdefault("min_cells_per_gene", 3)
    demo_config["qc"] = qc_cfg

    logging.info("Prepared synthetic datasets for demo (%d dataset(s))", len(dataset_entries))

    for dataset in dataset_entries:
        dataset_id = dataset["id"]
        with timer(f"QC ({dataset_id})"):
            qc.process_dataset_qc(dataset_id, demo_config)
        with timer(f"Doublets ({dataset_id})"):
            doublets.process_dataset_doublets(dataset_id, demo_config)

    integration_executed = False
    try:
        with timer("Integration"):
            integration.run_integration(demo_config)
            integration_executed = True
    except ImportError as exc:
        logging.warning("Integration dependency missing (%s); using fallback integration.", exc)
    except Exception as exc:  # pragma: no cover - defensive
        logging.warning("Integration failed (%s); using fallback integration.", exc)

    if not integration_executed:
        _run_fallback_integration(demo_config, [entry["id"] for entry in dataset_entries])

    with timer("Annotation"):
        annotate.run_annotation(demo_config)

    with timer("Export"):
        export.run_export(demo_config)

    with timer("Visualization"):
        viz.generate_all_figures(demo_config)

    demo_config_path = _write_demo_config_file(demo_config)
    with timer("Report"):
        generate_report(str(demo_config_path))

    logging.info("Demo pipeline completed. Artifacts are available under %s", Path("processed").absolute())


def create_minimal_demo(config_path: str) -> None:
    """Create lightweight demo assets without optional dependencies."""
    import anndata as ad
    import matplotlib.pyplot as plt

    setup_logging()
    config = load_config(config_path)
    set_seed(config.get("seed", 0))

    logging.info("Creating minimal demo assets...")

    ensure_dir(Path("processed"))
    ensure_dir(Path("processed/figures"))
    ensure_dir(Path("processed/cellxgene_release"))

    rng = np.random.default_rng(0)
    n_cells = 80
    n_genes = 60

    counts = rng.poisson(lam=2.0, size=(n_cells, n_genes)).astype(np.float32)
    counts += rng.binomial(1, 0.4, size=(n_cells, n_genes))

    obs = pd.DataFrame({
        "dataset_id": rng.choice([d["id"] for d in config.get("datasets", [])[:2]] or ["DEMO"], size=n_cells),
        "cancer_type": rng.choice(["Melanoma", "NSCLC", "Breast"], size=n_cells),
        "cell_type": rng.choice(["CD8_T", "CD4_T", "NK", "B_cell"], size=n_cells),
    })
    var = pd.DataFrame(index=[f"GENE_{i}" for i in range(n_genes)])

    adata = ad.AnnData(counts, obs=obs, var=var)
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names_make_unique()

    integrated_path = Path("processed") / "integrated_atlas.h5ad"
    annotated_path = Path("processed") / "integrated_annotated.h5ad"
    cellxgene_path = Path("processed") / "cellxgene_release" / "atlas.h5ad"

    adata.write_h5ad(integrated_path)

    annotated = adata.copy()
    annotated.write_h5ad(annotated_path)

    ensure_dir(cellxgene_path.parent)
    annotated.write_h5ad(cellxgene_path)

    plt.figure(figsize=(6, 4))
    obs["cell_type"].value_counts().plot(kind="bar", color="#4c72b0")
    plt.ylabel("Cells")
    plt.title("Cell Type Distribution (Demo)")
    plt.tight_layout()
    plt.savefig(Path("processed/figures") / "cell_type_distribution.png", dpi=300)
    plt.close()

    plt.figure(figsize=(6, 4))
    obs["cancer_type"].value_counts().plot(kind="bar", color="#55a868")
    plt.ylabel("Cells")
    plt.title("Cancer Type Distribution (Demo)")
    plt.tight_layout()
    plt.savefig(Path("processed/figures") / "cancer_type_distribution.png", dpi=300)
    plt.close()

    generate_report(config_path)

    logging.info("Minimal demo artifacts generated successfully.")


def _generate_tcr_section(config: Dict, metrics_dir: Path, figures_dir: Path, embed_figures: bool) -> List[str]:
    """Generate TCR/BCR repertoire section for report.

    Returns a list of markdown lines. Returns empty list if TCR is disabled or no data available.
    """
    import json
    import os

    lines: List[str] = []

    # Check if TCR is enabled
    receptor_cfg = config.get("tcr") or config.get("receptor", {}) or {}
    if not receptor_cfg.get("enabled", False):
        return lines

    # Determine TCR metrics directory
    tcr_metrics_dir = Path(
        receptor_cfg.get(
            "qc_metrics_dir",
            receptor_cfg.get("metrics_dir", "processed/metrics/tcr"),
        )
    )

    tcr_figures_dir = Path(
        receptor_cfg.get("figures_dir", "processed/figures/tcr")
    )

    # Check if any TCR data exists
    tcr_summary_path = tcr_metrics_dir / "tcr_summary.json"
    if not tcr_summary_path.exists():
        # TCR analysis hasn't run yet
        return lines

    # Start TCR section
    lines.append("## TCR/BCR Repertoire Analysis")
    lines.append("")

    try:
        tcr_summary = json.loads(tcr_summary_path.read_text())

        # Global summary metrics
        lines.append("### Global Repertoire Statistics")
        lines.append("")

        global_metrics = tcr_summary.get("global", {})
        if global_metrics:
            lines.append("| Metric | Value |")
            lines.append("|---|---|")

            if "total_cells" in global_metrics:
                lines.append(f"| Total cells with receptor data | {global_metrics['total_cells']:,} |")
            if "total_clonotypes" in global_metrics:
                lines.append(f"| Unique clonotypes | {global_metrics['total_clonotypes']:,} |")
            if "n_public_clonotypes" in global_metrics:
                lines.append(f"| Public clonotypes (shared across datasets) | {global_metrics['n_public_clonotypes']:,} |")

            # Repertoire overlap
            overlap_summary = global_metrics.get("repertoire_overlap_summary", {})
            if overlap_summary:
                if "mean_jaccard" in overlap_summary:
                    lines.append(f"| Mean repertoire overlap (Jaccard) | {overlap_summary['mean_jaccard']:.3f} |")
                if "mean_morisita_horn" in overlap_summary:
                    lines.append(f"| Mean repertoire similarity (Morisita-Horn) | {overlap_summary['mean_morisita_horn']:.3f} |")

            # CDR3 properties
            cdr3_props = global_metrics.get("cdr3_properties", {})
            if cdr3_props:
                if "mean_length" in cdr3_props:
                    lines.append(f"| Mean CDR3 length (AA) | {cdr3_props['mean_length']:.1f} |")
                if "mean_charge" in cdr3_props:
                    lines.append(f"| Mean CDR3 charge | {cdr3_props['mean_charge']:.2f} |")
                if "mean_hydrophobicity" in cdr3_props:
                    lines.append(f"| Mean CDR3 hydrophobicity | {cdr3_props['mean_hydrophobicity']:.2f} |")

            lines.append("")

        # Per-dataset metrics
        lines.append("### Per-Dataset Repertoire Metrics")
        lines.append("")

        dataset_metrics = tcr_summary.get("datasets", {})
        if dataset_metrics:
            lines.append("| Dataset | Cells | Clonotypes | Shannon Entropy | D50 | Top Clonotype Size |")
            lines.append("|---|---|---|---|---|---|")

            for dataset_id, metrics in sorted(dataset_metrics.items()):
                n_cells = metrics.get("n_cells", 0)
                n_clonotypes = metrics.get("n_clonotypes", 0)

                diversity = metrics.get("diversity", {})
                shannon = diversity.get("shannon_entropy", "N/A")
                shannon_str = f"{shannon:.2f}" if isinstance(shannon, (int, float)) else shannon

                d50 = diversity.get("D50", "N/A")
                d50_str = f"{int(d50)}" if isinstance(d50, (int, float)) else d50

                top_clonotypes = metrics.get("top_clonotypes", [])
                top_size = top_clonotypes[0]["cells"] if top_clonotypes else 0

                lines.append(
                    f"| {dataset_id} | {n_cells:,} | {n_clonotypes:,} | "
                    f"{shannon_str} | {d50_str} | {top_size:,} |"
                )

            lines.append("")

        # V/J gene usage summary
        lines.append("### V/J Gene Usage")
        lines.append("")
        lines.append("Top 10 most frequent V genes:")
        lines.append("")

        v_usage = global_metrics.get("v_gene_usage", {})
        if v_usage:
            lines.append("| V Gene | Cells |")
            lines.append("|---|---|")
            for gene, count in sorted(v_usage.items(), key=lambda x: x[1], reverse=True)[:10]:
                lines.append(f"| {gene} | {count:,} |")
            lines.append("")
        else:
            lines.append("_No V gene usage data available._")
            lines.append("")

        lines.append("Top 10 most frequent J genes:")
        lines.append("")

        j_usage = global_metrics.get("j_gene_usage", {})
        if j_usage:
            lines.append("| J Gene | Cells |")
            lines.append("|---|---|")
            for gene, count in sorted(j_usage.items(), key=lambda x: x[1], reverse=True)[:10]:
                lines.append(f"| {gene} | {count:,} |")
            lines.append("")
        else:
            lines.append("_No J gene usage data available._")
            lines.append("")

        # Public clonotypes
        public_clonotypes_path = tcr_metrics_dir / "public_clonotypes.json"
        if public_clonotypes_path.exists():
            public_data = json.loads(public_clonotypes_path.read_text())

            n_public = public_data.get("n_public_clonotypes", 0)
            if n_public > 0:
                lines.append("### Public Clonotypes")
                lines.append("")
                lines.append(f"Identified **{n_public:,}** clonotypes shared across multiple datasets.")
                lines.append("")

                top_public = public_data.get("top_public_clonotypes", [])[:5]
                if top_public:
                    lines.append("Top 5 most prevalent public clonotypes:")
                    lines.append("")
                    lines.append("| Clonotype ID | Datasets | Total Cells |")
                    lines.append("|---|---|---|")

                    for clono_info in top_public:
                        clono_id = clono_info["clonotype_id"]
                        n_datasets = clono_info["n_datasets"]
                        total_cells = clono_info["total_cells"]
                        datasets_str = ", ".join(clono_info["datasets"][:3])
                        if n_datasets > 3:
                            datasets_str += f" (+{n_datasets - 3} more)"

                        lines.append(f"| `{clono_id}` | {datasets_str} | {total_cells:,} |")

                    lines.append("")

        # Embed TCR figures
        if embed_figures:
            lines.append("### Repertoire Visualizations")
            lines.append("")

            tcr_figure_names = [
                ("clonotype_frequency_top20.png", "Top 20 Expanded Clonotypes"),
                ("repertoire_diversity_by_cancer_type.png", "Repertoire Diversity by Cancer Type"),
                ("umap_clonotype_expansion.png", "UMAP Colored by Clonotype Size"),
                ("cdr3_spectratype_by_chain.png", "CDR3 Length Distribution (Spectratype)"),
                ("repertoire_overlap_jaccard.png", "Repertoire Overlap (Jaccard Index)"),
                ("vj_pairing_heatmap.png", "V-J Gene Pairing Frequencies"),
            ]

            report_dir = Path("processed")
            figures_found = False

            for fig_name, title in tcr_figure_names:
                fig_path = tcr_figures_dir / fig_name
                if fig_path.exists():
                    rel_path = os.path.relpath(fig_path, report_dir)
                    lines.append(f"#### {title}")
                    lines.append(f"![{title}]({rel_path})")
                    lines.append("")
                    figures_found = True

            if not figures_found:
                lines.append("_TCR/BCR figures not yet generated._")
                lines.append("")

        # Additional metrics files
        lines.append("### Additional TCR/BCR Metrics")
        lines.append("")
        lines.append(f"- **Summary metrics**: `{tcr_metrics_dir}/tcr_summary.json`")
        lines.append(f"- **Repertoire overlap**: `{tcr_metrics_dir}/repertoire_overlap.json`")
        lines.append(f"- **Public clonotypes**: `{tcr_metrics_dir}/public_clonotypes.json`")

        # List per-dataset metrics
        dataset_metric_files = sorted(tcr_metrics_dir.glob("*_tcr_metrics.json"))
        if dataset_metric_files:
            lines.append(f"- **Per-dataset metrics**: {len(dataset_metric_files)} files in `{tcr_metrics_dir}/`")

        lines.append("")

    except Exception as e:
        # Gracefully degrade on error
        lines.append("_TCR/BCR analysis metrics could not be loaded._")
        lines.append(f"_Error: {str(e)}_")
        lines.append("")

    return lines


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

    # Generate TCR/BCR section (replaces old fragment-based approach)
    tcr_section = _generate_tcr_section(config, metrics_dir, figures_dir, embed_figures)
    if tcr_section:
        lines.extend(tcr_section)

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
