import yaml
from pathlib import Path

from src.atlas.workflow import collect_pipeline_targets, dataset_ids

# Load configuration
with open("config/atlas.yaml", "r") as f:
    config = yaml.safe_load(f)

# Extract dataset IDs
DATASETS = dataset_ids(config)

PIPELINE_TARGETS = collect_pipeline_targets(config)

METRICS_DIR = Path(config["outputs"].get("metrics_dir", "processed/metrics"))
FIGURES_DIR = Path(config["outputs"].get("figures_dir", "processed/figures"))
QC_PLOTS_DIR = Path(config.get("qc", {}).get("qc_plots_dir", FIGURES_DIR / "qc"))
DOUBLETS_METRICS_DIR = METRICS_DIR / "doublets"
DOUBLETS_FIGURES_DIR = FIGURES_DIR / "doublets"

rule all:
    input:
        PIPELINE_TARGETS

rule download_raw:
    output:
        expand("data/raw/{dataset}_downloaded.flag", dataset=DATASETS)
    params:
        config_file="config/atlas.yaml"
    shell:
        "python -m src.atlas.io --download --config {params.config_file}"

rule qc_one:
    input:
        "data/raw/{dataset}_downloaded.flag"
    output:
        filtered="data/interim/{dataset}.h5ad",
        qc_table="data/interim/{dataset}_qc.tsv",
        qc_summary=lambda wildcards: str(METRICS_DIR / f"{wildcards.dataset}_qc_summary.json"),
        qc_plot_umi=lambda wildcards: str(QC_PLOTS_DIR / f"{wildcards.dataset}_umi_violin.png"),
        qc_plot_mt=lambda wildcards: str(QC_PLOTS_DIR / f"{wildcards.dataset}_mt_violin.png")
    params:
        config_file="config/atlas.yaml",
        dataset="{dataset}"
    shell:
        "python -m src.atlas.qc --dataset {params.dataset} --config {params.config_file}"

rule doublets_one:
    input:
        "data/interim/{dataset}.h5ad"
    output:
        filtered="data/interim/{dataset}.doublet_filtered.h5ad",
        annotated="data/interim/{dataset}.with_doublets.h5ad",
        summary=lambda wildcards: str(DOUBLETS_METRICS_DIR / f"{wildcards.dataset}_doublet_summary.json"),
        histogram=lambda wildcards: str(DOUBLETS_FIGURES_DIR / f"{wildcards.dataset}_scrublet_hist.png")
    params:
        config_file="config/atlas.yaml",
        dataset="{dataset}"
    shell:
        "python -m src.atlas.doublets --dataset {params.dataset} --config {params.config_file}"

rule integrate:
    input:
        expand("data/interim/{dataset}.doublet_filtered.h5ad", dataset=DATASETS)
    output:
        atlas="processed/integrated_atlas.h5ad",
        metrics="processed/metrics/integration_metrics.json"
    params:
        config_file="config/atlas.yaml"
    shell:
        "python -m src.atlas.integration --config {params.config_file}"

rule annotate:
    input:
        "processed/integrated_atlas.h5ad"
    output:
        atlas="processed/integrated_annotated.h5ad",
        summary="processed/metrics/annotation_summary.json",
        scores="processed/metrics/annotation_scores.tsv",
        confusion="processed/metrics/annotation_confusion.tsv"
    params:
        config_file="config/atlas.yaml"
    shell:
        "python -m src.atlas.annotate --config {params.config_file}"

rule export_cellxgene:
    input:
        "processed/integrated_annotated.h5ad"
    output:
        "processed/cellxgene_release/atlas.h5ad"
    params:
        config_file="config/atlas.yaml"
    shell:
        "python -m src.atlas.export --cellxgene --config {params.config_file}"

rule figures:
    input:
        "processed/integrated_annotated.h5ad"
    output:
        "processed/figures/umap_by_cell_type.png",
        "processed/figures/umap_by_dataset.png",
        "processed/figures/umap_by_cancer_type.png",
        "processed/figures/proportions_by_cancer_type.png"
    params:
        config_file="config/atlas.yaml"
    shell:
        "python -m src.atlas.viz --all --config {params.config_file}"

rule benchmark:
    input:
        "processed/integrated_annotated.h5ad"
    output:
        str(METRICS_DIR / "benchmarking.json")
    params:
        config_file="config/atlas.yaml"
    run:
        if not config.get("benchmarking", {}).get("enabled", False):
            open(output[0], "w").write("{}\n")
        else:
            shell("python -m src.atlas.benchmark --config {params.config_file}")

rule report:
    input:
        "processed/integrated_annotated.h5ad",
        "processed/figures/umap_by_cell_type.png"
    output:
        "processed/report.md"
    params:
        config_file="config/atlas.yaml"
    shell:
        "python -c \"from src.atlas.utils import generate_report; generate_report('{params.config_file}')\""
