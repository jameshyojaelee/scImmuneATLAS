# scImmuneATLAS

[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

A comprehensive, production-ready single-cell immune atlas for tumor-infiltrating lymphocytes (TILs).

scImmuneATLAS integrates multiple public scRNA-seq datasets, providing a standardized, reproducible pipeline for immune cell analysis in cancer contexts.

## ✨ Key Features

- Dataset ingest + fetch: local H5AD support and optional CELLxGENE Census downloader
- **Data integrity guarantees**: Pandera schema validation for AnnData objects and optional SHA256 checksum verification for downloads
- QC dashboards: canonical gene/count thresholds, per-dataset summaries, violin plots saved to `processed/figures/qc/`
- Doublet diagnostics: single-pass Scrublet with histograms + JSON summaries in `processed/metrics/doublets/`
- Batch correction: shared-HVG integration with Harmony or scVI (GPU auto-detect) and mixing metrics (`processed/metrics/integration_metrics.json`)
- Annotation: direction-aware marker scoring, confidence gating, optional reference fallback, and confusion matrices
- Reporting: Markdown report embedding QC/doublet/integration/annotation metrics and key figures
- Tooling: Snakemake pipeline, typed CLI (`scimmuneatlas`), unit tests, and portfolio-ready notebooks/docs

## 🔄 Analysis Workflow

1. Data ingestion
   - H5AD directly; MTX+genes/barcodes supported
   - Optional fetch from CELLxGENE Census
2. QC + doublet removal
   - Applies gene/count/mito thresholds, exports per-dataset QC tables/plots (`data/interim/*_qc.tsv`, `processed/figures/qc/`)
   - Runs Scrublet once per dataset, writes annotated + filtered H5ADs and metrics (`processed/metrics/doublets/`)
3. Batch integration
   - Harmony or scVI (shared HVGs, GPU auto-detect)
   - Output: `processed/integrated_atlas.h5ad` + integration diagnostics (`processed/metrics/integration_metrics.json`)
4. Cell type annotation
   - Marker-direction-aware scoring with confidence gating and optional reference fallback
   - Output: `processed/integrated_annotated.h5ad`, annotation scores, confusion matrices
5. Visualization, benchmarking & export
   - Figures, Markdown report with embedded metrics, optional benchmarking against reference atlas, cellxgene export

## 🚀 Quick Start

### Prerequisites
- Python 3.10+
- Conda/Mamba (or micromamba)

### Installation

1) Clone
```
git clone https://github.com/YOUR_ORG_OR_USERNAME/scImmuneATLAS.git
cd scImmuneATLAS
```

2) Create environment (choose one)

- **Conda / Mamba**
  ```
  mamba env create -f env.yml   # GPU-enabled if CUDA is available
  mamba activate immune-atlas
  # CPU-only alternative:
  # conda env create -f env.cpu.yml && conda activate immune-atlas
  # No micromamba? use env.nomamba.yml
  ```
- **Python venv**
  ```
  python -m venv .venv
  source .venv/bin/activate
  pip install -e .[dev] scanpy scrublet harmonypy scvi-tools scib matplotlib seaborn
  ```

3) Dev tooling
```
pre-commit install
```

### Run Demo

End-to-end synthetic demo:
```
make demo
```
Outputs go to `processed/figures/` and `processed/cellxgene_release/`.

## Data Acquisition

Use the example `.h5ad` inputs in `data/raw/` referenced by `config/atlas.yaml`, or fetch new subsets from CELLxGENE Census.

- Provided inputs: ensure paths in `config/atlas.yaml` point to your local files in `data/raw/`.
- Fetch from CELLxGENE Census (network required):
```
make env-bootstrap             # optional: bootstrap micromamba locally
make fetch-no-cap              # runs scripts/fetch_cellxgene_no_cap.sh (parallel jobs)
make run-after-fetch           # or: bash scripts/run_pipeline_after_fetch.sh
# Alternatively: scripts/watch_fetch_then_pipeline.sh to auto-start after fetch
```

## Datasets Included

Configured datasets (see `config/atlas.yaml`). The first three are present locally; the others typically require fetching or adding matching `.h5ad` files under `data/raw/`.

- MELANOMA_CENSUS: Melanoma tumor microenvironment immune scRNA-seq subset  
  File: `data/raw/melanoma_census.h5ad`
- NSCLC_CENSUS: Non-small cell lung carcinoma (NSCLC) immune scRNA-seq subset  
  File: `data/raw/nsclc_census.h5ad`
- BREAST_CENSUS: Breast carcinoma immune scRNA-seq subset  
  File: `data/raw/breast_census.h5ad`
- NSCLC_PEMBROLIZUMAB_CENSUS: NSCLC cohort with pembrolizumab (anti–PD-1) therapy  
  File: `data/raw/nsclc_pembro_census.h5ad`
- NSCLC_NIVO_IPI_CENSUS: NSCLC cohort treated with nivolumab + ipilimumab (PD-1 + CTLA-4)  
  File: `data/raw/nsclc_nivo_ipi_census.h5ad`
- RCC_COMBINATION_CENSUS: Renal cell carcinoma combination therapy cohorts (immune checkpoint inhibitors)  
  File: `data/raw/rcc_combo_census.h5ad`

## Project Structure

```
scImmuneATLAS/
├── README.md                    # This file
├── LICENSE                      # MIT license
├── env.yml                      # Main Conda/mamba environment (GPU optional)
├── env.cpu.yml                  # CPU-only environment
├── env.nomamba.yml              # Fallback env without micromamba helpers
├── Makefile                     # Lint, test, demo, app, fetch, etc.
├── Snakefile                    # Snakemake workflow
├── config/
│   └── atlas.yaml               # Analysis configuration
├── data/
│   ├── raw/                     # Original datasets (inputs)
│   ├── interim/                 # Intermediate processing
│   └── external/                # Reference metadata (e.g., metadata.tsv)
├── processed/
│   ├── figures/                 # Publication-ready figures (UMAPs, QC, doublets)
│   ├── metrics/                 # JSON/TSV diagnostics (QC, doublets, integration, annotation, benchmarking)
│   └── cellxgene_release/       # Viewer-compatible outputs
├── notebooks/                   # Interpretation + metadata notebooks
├── src/atlas/                   # Core analysis modules
│   ├── io.py, qc.py, doublets.py, integration.py, annotate.py, benchmark.py, viz.py, export.py, cli.py
│   └── markers/immune_markers_human.tsv
├── app/                         # Streamlit web app
│   └── streamlit_app.py
├── scripts/                     # Helper scripts
│   ├── fetch_cellxgene_no_cap.sh
│   ├── run_pipeline_after_fetch.sh
│   └── watch_fetch_then_pipeline.sh
├── tests/                       # Unit tests
└── logs/                        # Logs from fetch/pipeline
```

## ⚙️ Configuration

Edit `config/atlas.yaml` to customize analysis. Example:
```
datasets:
  - id: "NSCLC_CENSUS"
    url: "data/raw/nsclc_census.h5ad"
    cancer_type: "NSCLC"
  - id: "MELANOMA_CENSUS"
    url: "data/raw/melanoma_census.h5ad"
    cancer_type: "Melanoma"

integration:
  method: "scvi"           # or "harmony"
  batch_key: "dataset_id"

annotation:
  marker_genes: "src/atlas/markers/immune_markers_human.tsv"
```

## 🏃‍♂️ Running the Pipeline

### Full pipeline
```
snakemake -j 8          # recommended
make all                # Make wrapper
```

The workflow materializes diagnostics throughout (`processed/metrics/`, `processed/figures/qc/`, `processed/metrics/doublets/`) and fails fast if expected artifacts are missing.

### CLI shortcuts
Install the package (editable mode) and drive each stage via the CLI:
```
scimmuneatlas validate-data --config config/atlas.yaml   # Validate dataset integrity
scimmuneatlas qc --config config/atlas.yaml
scimmuneatlas doublets --config config/atlas.yaml
scimmuneatlas integrate --config config/atlas.yaml
scimmuneatlas annotate --config config/atlas.yaml
scimmuneatlas viz --config config/atlas.yaml --log-level INFO
scimmuneatlas export --config config/atlas.yaml --cellxgene
scimmuneatlas benchmark --config config/atlas.yaml
scimmuneatlas report --config config/atlas.yaml
```

### Individual modules (optional)
```
# Download (write metadata flags and a dataset table)
python -m src.atlas.io --download --config config/atlas.yaml

# QC (all datasets) and doublets
python -m src.atlas.qc --run --config config/atlas.yaml
python -m src.atlas.doublets --run --config config/atlas.yaml

# Integration (override method if desired)
python -m src.atlas.integration --method harmony --config config/atlas.yaml

# Annotation, export, figures, benchmarking, report
python -m src.atlas.annotate --config config/atlas.yaml
python -m src.atlas.viz --all --config config/atlas.yaml
python -m src.atlas.export --cellxgene --config config/atlas.yaml
python -m src.atlas.benchmark --config config/atlas.yaml
python -c "from src.atlas.utils import generate_report; generate_report('config/atlas.yaml')"
```

## 📤 Outputs

- `processed/integrated_atlas.h5ad`: Integrated atlas after batch correction
- `processed/integrated_annotated.h5ad`: Annotated atlas with immune cell types
- `processed/cellxgene_release/atlas.h5ad`: cellxgene-compatible release
- `processed/figures/`: Publication-ready visualizations (UMAPs, QC violins, Scrublet histograms)
- `processed/metrics/`: QC summaries, doublet metrics, integration diagnostics, annotation scores, benchmarking results
- `processed/report.md`: Analysis summary report with embedded tables and figures

## 🖥️ Interactive Explorer

Launch Streamlit app:
```
make app
```

Features: UMAP exploration, filtering, gene expression overlays, proportions, dataset summary.

## 🧪 Development

- Lint: `make lint`
- Tests: `make test`
- Validate data: `make validate-data`
- Pre-commit hooks: `pre-commit install`
- Local test run with venv: `source .venv/bin/activate && PYTHONPATH=$PWD/src pytest -q`

### Data Validation & Integrity

The pipeline includes built-in data validation to ensure dataset quality:

**Schema Validation**: All datasets are validated against a Pandera schema that checks:
- Required metadata columns (`dataset_id`, `cancer_type`, `platform`)
- Non-empty observations and variables
- Unique cell and gene identifiers

**Checksum Verification**: Optional SHA256 hash verification for downloaded files. Add hashes to `config/atlas.yaml`:

```yaml
datasets:
  - id: "EXAMPLE_DATASET"
    url: "https://example.com/data.h5ad"
    url_sha256: "abc123def456..."  # Optional SHA256 hash
    cancer_type: "melanoma"
    platform: "10x"
```

**Validation CLI**: Run `make validate-data` or `scimmuneatlas validate-data` to check all datasets. The command:
- Verifies file existence
- Loads each dataset
- Validates schema compliance
- Reports actionable error messages
- Exits with non-zero code on failures (CI-friendly)

## ❓ FAQ & Troubleshooting

- Memory: use backed mode for large datasets
  ```python
  import anndata as ad
  adata = ad.read_h5ad("large_file.h5ad", backed='r')
  ```
- GPU: scVI uses GPU if available
  ```
  mamba install pytorch-cuda -c pytorch -c nvidia
  ```
- Env isolation: user site-packages are disabled to avoid mixing with `~/.local`; prefer running within the provided Conda/mamba environment.
- Integration choice: Harmony is lighter/faster; scVI handles complex batch effects better. Both stages share HVGs so diagnostics in `processed/metrics/integration_metrics.json` are comparable.
- Missing metrics? Run stages through Snakemake or the CLI so QC/doublet/integration/annotation summaries populate `processed/metrics/`.

## 📓 Portfolio Assets

- `notebooks/atlas_interpretation.ipynb`: Guided notebook for showcasing biological insights
- `notebooks/metadata_integration_example.ipynb`: Template for merging clinical/TCR metadata
- `docs/video_walkthrough_script.md`: Script outline for a short project walkthrough
- `docs/blog_post_outline.md`: Blog post scaffold detailing design decisions
- `docs/portfolio_one_pager.md`: Printable summary for interviews and lab meetings

## 📚 Citation

```
@software{scImmuneATLAS,
  title = {scImmuneATLAS: A comprehensive single-cell immune atlas framework},
  author = {Your Name},
  year = {2024},
  url = {https://github.com/YOUR_ORG_OR_USERNAME/scImmuneATLAS}
}
```

## 🤝 Contributing

1) Fork the repository  
2) Create a feature branch (`git checkout -b feature/amazing-feature`)  
3) Run tests (`make test`) and linting (`make lint`)  
4) Commit (`git commit -m 'Add amazing feature'`)  
5) Push (`git push origin feature/amazing-feature`)  
6) Open a Pull Request

## 📄 License

MIT License — see `LICENSE`.
