import yaml
from pathlib import Path

# Load configuration
with open("config/atlas.yaml", "r") as f:
    config = yaml.safe_load(f)

# Extract dataset IDs
DATASETS = [d["id"] for d in config["datasets"]]

rule all:
    input:
        "processed/integrated_annotated.h5ad",
        "processed/cellxgene_release/atlas.h5ad",
        "processed/figures/umap_by_cell_type.png",
        "processed/figures/umap_by_dataset.png",
        "processed/figures/umap_by_cancer_type.png",
        "processed/figures/proportions_by_cancer_type.png",
        "processed/report.md"

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
        "data/interim/{dataset}.h5ad"
    params:
        config_file="config/atlas.yaml",
        dataset="{dataset}"
    shell:
        "python -m src.atlas.qc --dataset {params.dataset} --config {params.config_file}"

rule doublets_one:
    input:
        "data/interim/{dataset}.h5ad"
    output:
        "data/interim/{dataset}.doublet_filtered.h5ad"
    params:
        config_file="config/atlas.yaml",
        dataset="{dataset}"
    shell:
        "python -m src.atlas.doublets --dataset {params.dataset} --config {params.config_file}"

rule integrate:
    input:
        expand("data/interim/{dataset}.doublet_filtered.h5ad", dataset=DATASETS)
    output:
        "processed/integrated_atlas.h5ad"
    params:
        config_file="config/atlas.yaml"
    shell:
        "python -m src.atlas.integration --config {params.config_file}"

rule annotate:
    input:
        "processed/integrated_atlas.h5ad"
    output:
        "processed/integrated_annotated.h5ad"
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
