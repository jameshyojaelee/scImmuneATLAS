import yaml
from pathlib import Path

from src.atlas.workflow import collect_pipeline_targets, dataset_ids
from src.atlas import receptor as receptor_module

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
RECEPTOR_ENABLED = receptor_module.is_enabled(config)
RECEPTOR_DATASETS = list(receptor_module.datasets_with_receptor(config)) if RECEPTOR_ENABLED else []
_TCR_CFG = config.get("tcr") or config.get("receptor", {}) or {}
RECEPTOR_METRICS_DIR = Path(
    _TCR_CFG.get(
        "qc_metrics_dir",
        _TCR_CFG.get("metrics_dir", "processed/metrics/tcr"),
    )
)
RECEPTOR_FIGURES_DIR = Path(
    _TCR_CFG.get("figures_dir", "processed/figures/tcr")
)


def _dataset_config(dataset_id):
    for entry in config.get("datasets", []):
        if entry["id"] == dataset_id:
            return entry
    raise KeyError(f"Unknown dataset: {dataset_id}")


RECEPTOR_PARQUET_MAP = (
    {
        dataset: str(
            receptor_module.ingest.resolved_output_path(_dataset_config(dataset), config)
        )
        for dataset in RECEPTOR_DATASETS
    }
    if RECEPTOR_DATASETS
    else {}
)

RECEPTOR_QC_MAP = (
    {dataset: str(receptor_module.qc.get_output_path(dataset, config)) for dataset in RECEPTOR_DATASETS}
    if RECEPTOR_DATASETS
    else {}
)

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

if RECEPTOR_ENABLED and RECEPTOR_DATASETS:
    rule receptor_ingest:
        input:
            "data/raw/{dataset}_downloaded.flag"
        output:
            parquet=lambda wildcards: RECEPTOR_PARQUET_MAP[wildcards.dataset]
        params:
            config_file="config/atlas.yaml",
            dataset="{dataset}"
        shell:
            "python -m src.atlas.receptor --stage ingest --dataset {params.dataset} --config {params.config_file}"

    rule receptor_qc:
        input:
            parquet=lambda wildcards: RECEPTOR_PARQUET_MAP[wildcards.dataset]
        output:
            qc=lambda wildcards: RECEPTOR_QC_MAP[wildcards.dataset]
        params:
            config_file="config/atlas.yaml",
            dataset="{dataset}"
        shell:
            "python -m src.atlas.receptor --stage qc --dataset {params.dataset} --config {params.config_file}"

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

if RECEPTOR_ENABLED and RECEPTOR_DATASETS:
    rule receptor_analytics:
        input:
            atlas="processed/integrated_annotated.h5ad",
            qc=[RECEPTOR_QC_MAP[d] for d in RECEPTOR_DATASETS],
            tables=[RECEPTOR_PARQUET_MAP[d] for d in RECEPTOR_DATASETS]
        output:
            summary=str(RECEPTOR_METRICS_DIR / "repertoire_summary.json"),
            expansion=str(RECEPTOR_METRICS_DIR / "clonotype_expansion.tsv"),
            diversity=str(RECEPTOR_METRICS_DIR / "diversity_indices.json"),
            vj=str(RECEPTOR_METRICS_DIR / "vj_usage.tsv"),
            pairing=str(RECEPTOR_METRICS_DIR / "chain_pairing_summary.json"),
            v_usage_json=str(RECEPTOR_METRICS_DIR / "v_usage_summary.json"),
            report_fragment=str(RECEPTOR_METRICS_DIR / "report_section.md"),
            freq_fig=str(RECEPTOR_FIGURES_DIR / "clonotype_frequency.png"),
            umap_fig=str(RECEPTOR_FIGURES_DIR / "umap_clonal_expansion.png"),
            spectratype_fig=str(RECEPTOR_FIGURES_DIR / "cdr3_spectratype.png"),
            vj_fig=str(RECEPTOR_FIGURES_DIR / "vj_usage_heatmap.png"),
        params:
            config_file="config/atlas.yaml"
        shell:
            "python -m src.atlas.receptor --stage analytics --config {params.config_file}"

    rule tcr:
        input:
            atlas="processed/integrated_annotated.h5ad",
            doublet_filtered=expand("data/interim/{dataset}.doublet_filtered.h5ad", dataset=RECEPTOR_DATASETS),
            tables=[RECEPTOR_PARQUET_MAP[d] for d in RECEPTOR_DATASETS]
        output:
            summary=str(RECEPTOR_METRICS_DIR / "tcr_summary.json"),
            overlap=str(RECEPTOR_METRICS_DIR / "repertoire_overlap.json"),
            public=str(RECEPTOR_METRICS_DIR / "public_clonotypes.json"),
            # Main visualization outputs
            freq_top20=str(RECEPTOR_FIGURES_DIR / "clonotype_frequency_top20.png"),
            diversity_fig=str(RECEPTOR_FIGURES_DIR / "repertoire_diversity_by_cancer_type.png"),
            umap_expansion=str(RECEPTOR_FIGURES_DIR / "umap_clonotype_expansion.png"),
            spectratype_chain=str(RECEPTOR_FIGURES_DIR / "cdr3_spectratype_by_chain.png"),
            overlap_jaccard=str(RECEPTOR_FIGURES_DIR / "repertoire_overlap_jaccard.png"),
            vj_heatmap=str(RECEPTOR_FIGURES_DIR / "vj_pairing_heatmap.png"),
        params:
            config_file="config/atlas.yaml"
        shell:
            "python -m src.atlas.tcr --config {params.config_file}"

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

REPORT_INPUTS = [
    "processed/integrated_annotated.h5ad",
    "processed/figures/umap_by_cell_type.png",
]

if RECEPTOR_ENABLED and RECEPTOR_DATASETS:
    REPORT_INPUTS.append(str(RECEPTOR_METRICS_DIR / "repertoire_summary.json"))
    REPORT_INPUTS.append(str(RECEPTOR_METRICS_DIR / "report_section.md"))
    REPORT_INPUTS.append(str(RECEPTOR_METRICS_DIR / "tcr_summary.json"))


rule report:
    input:
        REPORT_INPUTS
    output:
        "processed/report.md"
    params:
        config_file="config/atlas.yaml"
    shell:
        "python -c \"from src.atlas.utils import generate_report; generate_report('{params.config_file}')\""
