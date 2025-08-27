# scImmuneATLAS

[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

**A comprehensive, production-ready single-cell immune atlas for tumor-infiltrating lymphocytes (TILs)**

scImmuneATLAS is a robust computational framework that integrates multiple public scRNA-seq datasets of tumor-infiltrating lymphocytes, providing researchers with a standardized, reproducible pipeline for immune cell analysis in cancer contexts.

## ✨ Key Features

- 🔬 **Multi-dataset Integration**: Seamlessly combine diverse TIL scRNA-seq datasets
- 🎯 **Advanced QC Pipeline**: Automated quality control and doublet detection
- 🧩 **Batch Effect Correction**: Support for both Harmony and scVI integration methods  
- 🏷️ **Immune Cell Annotation**: Marker-based cell type identification
- 📊 **Interactive Visualization**: Built-in Streamlit app and cellxgene compatibility
- 🔄 **Reproducible Workflow**: Snakemake-based pipeline with comprehensive testing

## 🔄 Analysis Workflow

The scImmuneATLAS pipeline consists of the following main steps:

1. **Raw Data Import**  
   - Accepts MTX, H5AD, or TSV files

2. **Quality Control & Doublet Removal**  
   - Filters low-quality cells and detects doublets  
   - Output: `data/interim/`

3. **Batch Integration**  
   - Integrates multiple datasets using Harmony or scVI  
   - Output: `processed/integrated_atlas.h5ad`

4. **Cell Type Annotation**  
   - Assigns immune cell types using marker-based scoring  
   - Output: `processed/integrated_annotated.h5ad`

5. **Visualization & Export**  
   - Generates figures, reports, and interactive viewers  
   - Output: Interactive dashboard, figures, and reports

**Workflow Overview:**

## 📊 Output Gallery

scImmuneATLAS generates publication-ready visualizations including:

- **🗺️ UMAP Embeddings**: Interactive plots colored by cell type, dataset, and cancer type
- **📈 Cell Composition**: Stacked barplots showing immune cell proportions across conditions  
- **🔥 Marker Heatmaps**: Expression profiles for immune cell type signatures
- **📱 Interactive Dashboard**: Streamlit-based explorer for real-time data interaction

## 🚀 Quick Start

### Prerequisites
- Python 3.8+
- Conda/Mamba package manager

### Installation

1. **Clone the repository**
```bash
git clone https://github.com/YOUR_USERNAME/scImmuneATLAS.git
cd scImmuneATLAS
```

2. **Set up environment**
```bash
mamba env create -f env.yml
mamba activate immune-atlas
```

3. **Install development tools**
```bash
pre-commit install
```

### Run Demo

Test the complete pipeline with synthetic data:

```bash
make demo
```

This command will:
- ✅ Generate synthetic TIL data
- ✅ Perform quality control and doublet detection  
- ✅ Integrate datasets using batch correction
- ✅ Annotate immune cell types
- ✅ Generate visualizations and export results

Results will be saved in `processed/figures/` and `processed/cellxgene_release/`.

## 📁 Project Structure

```
scImmuneATLAS/
├── 📋 README.md              # This file
├── 📜 LICENSE                # MIT license
├── 🐍 env.yml               # Conda environment
├── ⚙️ Makefile              # Build automation
├── 🔧 Snakefile             # Workflow definition
├── 📊 config/
│   └── atlas.yaml           # Analysis parameters
├── 💾 data/                 # Data directory
│   ├── raw/                 # Original datasets
│   ├── interim/             # Intermediate processing
│   ├── processed/           # Final outputs
│   └── external/            # Reference data
├── 📓 notebooks/            # Analysis notebooks
│   ├── 00_explore_inputs.ipynb
│   ├── 01_qc_and_doublets.ipynb
│   ├── 02_integration.ipynb
│   ├── 03_annotation_and_markers.ipynb
│   └── 04_proportions_and_signatures.ipynb
├── 🔬 src/atlas/           # Core analysis modules
│   ├── io.py               # Data I/O operations
│   ├── qc.py               # Quality control
│   ├── integration.py      # Batch correction
│   ├── annotate.py         # Cell type annotation
│   ├── viz.py              # Visualization
│   └── markers/            # Reference gene sets
├── 🖥️ app/                 # Interactive applications
│   └── streamlit_app.py    # Web dashboard
└── 🧪 tests/               # Unit tests
    ├── test_smoke.py
    ├── test_io.py
    └── test_qc.py
```

## ⚙️ Configuration

Edit `config/atlas.yaml` to customize your analysis:

```yaml
# Example configuration
datasets:
  - name: "dataset1"
    url: "path/to/data.h5ad"
    cancer_type: "NSCLC"
  - name: "dataset2" 
    url: "https://example.com/data.h5ad"
    cancer_type: "Melanoma"

integration:
  method: "harmony"  # or "scvi"
  batch_key: "dataset"

annotation:
  marker_genes: "src/atlas/markers/immune_markers_human.tsv"
```

## 🏃‍♂️ Running the Pipeline

### Full Pipeline
```bash
# Using Snakemake (recommended)
snakemake -j 8

# Using Make wrapper  
make all
```

### Individual Steps
```bash
# Quality control
python -m src.atlas.qc --config config/atlas.yaml

# Integration  
python -m src.atlas.integration --method harmony --config config/atlas.yaml

# Annotation
python -m src.atlas.annotate --config config/atlas.yaml
```

## 📤 Outputs

| File | Description |
|------|-------------|
| `processed/integrated_atlas.h5ad` | 🧬 Integrated atlas after batch correction |
| `processed/integrated_annotated.h5ad` | 🏷️ Annotated atlas with immune cell types |
| `processed/cellxgene_release/atlas.h5ad` | 🌐 cellxgene-compatible release |
| `processed/figures/` | 📊 Publication-ready visualizations |
| `processed/report.md` | 📋 Analysis summary report |

## 🖥️ Interactive Explorer

Launch the web-based data explorer:

```bash
make app
```

Features:
- 🔍 Interactive UMAP exploration
- 📊 Real-time filtering and selection
- 📈 Gene expression visualization  
- 💾 Export selected data



## ❓ FAQ & Troubleshooting

<details>
<summary><strong>Memory Issues</strong></summary>

For large datasets, use backed mode to reduce memory usage:
```python
adata = anndata.read_h5ad("large_file.h5ad", backed='r')
```
</details>

<details>
<summary><strong>GPU Support</strong></summary>

scVI will automatically detect GPU if available. To enable:
```bash
mamba install pytorch-cuda -c pytorch -c nvidia
```
</details>

<details>
<summary><strong>Integration Methods</strong></summary>

- **Harmony**: Faster, lighter weight, good for most use cases
- **scVI**: Better for complex batch effects, requires more compute
</details>

<details>
<summary><strong>Reproducibility</strong></summary>

Random seeds are set for deterministic results. Minor UMAP variations may still occur due to numerical precision.
</details>

## 📚 Citation

If you use scImmuneATLAS in your research, please cite:

```bibtex
@software{scImmuneATLAS,
  title = {scImmuneATLAS: A comprehensive single-cell immune atlas framework},
  author = {Your Name},
  year = {2024},
  url = {https://github.com/YOUR_USERNAME/scImmuneATLAS}
}
```

## 🤝 Contributing

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Run tests (`make test`) and linting (`make lint`)
4. Commit your changes (`git commit -m 'Add amazing feature'`)
5. Push to the branch (`git push origin feature/amazing-feature`)
6. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

<p align="center">
  <strong>🧬 Built for the single-cell community by researchers, for researchers 🔬</strong>
</p>
