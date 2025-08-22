Single-cell Immune Atlas
========================

Project goal
------------
Build a production-quality, reproducible Single-cell Immune Atlas that integrates multiple public tumor-infiltrating lymphocyte (TIL) scRNA-seq datasets, performs QC and doublet removal, corrects batch effects (Harmony or scVI), annotates immune cell types via markers, and publishes an interactive viewer (cellxgene-ready h5ad and a minimal Streamlit browser).

Pipeline overview
-----------------

```
Raw data (MTX/H5AD/TSV)
   │
   ├── download_raw  → data/raw/
   │
   ├── qc_one        → data/interim/{dataset}.h5ad
   │
   ├── doublets_one  → data/interim/{dataset}.doublet_filtered.h5ad
   │
   ├── integrate (Harmony | scVI) → processed/integrated_atlas.h5ad
   │
   ├── annotate (marker-based)    → processed/integrated_annotated.h5ad
   │
   ├── export_cellxgene           → processed/cellxgene_release/atlas.h5ad
   │
   └── figures + report           → processed/figures/, processed/report.md
```

Example figures
---------------
This repository generates:
- UMAP colored by cell type, dataset, and cancer type
- Stacked barplots of cell-type proportions per cancer type
- A marker heatmap for a small panel of immune reference genes

Quickstart
----------

1) Create and activate environment

```bash
mamba env create -f env.yml
mamba activate immune-atlas
```

2) Install pre-commit hooks

```bash
pre-commit install
```

3) Run the tiny end-to-end demo

```bash
make demo
```

This will synthesize a small dataset, run QC, integration, marker-based annotation, export a cellxgene-compatible h5ad, and save figures under `processed/figures/`.

Repository layout
-----------------

```
single-cell-immune-atlas/
  README.md
  LICENSE
  env.yml
  Makefile
  .gitignore
  .pre-commit-config.yaml
  pyproject.toml
  setup.cfg
  Snakefile
  config/
    atlas.yaml
  data/
    raw/
    interim/
    processed/
    external/
  notebooks/
    00_explore_inputs.ipynb
    01_qc_and_doublets.ipynb
    02_integration.ipynb
    03_annotation_and_markers.ipynb
    04_proportions_and_signatures.ipynb
  src/atlas/
    __init__.py
    io.py
    qc.py
    doublets.py
    integration.py
    annotate.py
    viz.py
    export.py
    utils.py
    markers/
      immune_markers_human.tsv
  app/
    streamlit_app.py
  tests/
    test_smoke.py
    test_io.py
    test_qc.py
```

Configuration
-------------
Datasets and key parameters live in `config/atlas.yaml`. Replace the placeholder URLs with real URLs or local file paths when running on real data. The code accepts local paths or HTTP(S) links.

How to run
----------
- Snakemake pipeline:

```bash
snakemake -j 8
```

- Or use GNU Make as a wrapper:

```bash
make all
```

Outputs
-------
- `processed/integrated_atlas.h5ad` — integrated atlas after batch correction
- `processed/integrated_annotated.h5ad` — annotated atlas with immune cell types
- `processed/cellxgene_release/atlas.h5ad` — cellxgene-compatible release object
- `processed/figures/` — static UMAPs, barplots, and heatmaps
- `processed/report.md` — lightweight markdown summary of the run

Interactive viewer
------------------
Launch the minimal Streamlit app:

```bash
make app
```

The app reads `processed/integrated_annotated.h5ad` and the configuration in `config/atlas.yaml`.

Daily Plan
----------

| Day | Task | Command(s) | Output |
| --- | --- | --- | --- |
| Mon | Download datasets & write metadata table | `python -m src.atlas.io --download --config config/atlas.yaml` | `data/raw/*`, `data/external/metadata.tsv` |
| Tue | QC + doublets per dataset | `python -m src.atlas.qc --run --config config/atlas.yaml` then `python -m src.atlas.doublets --run --config config/atlas.yaml` | `data/interim/*.doublet_filtered.h5ad` |
| Wed | Integration (scVI or Harmony) | `python -m src.atlas.integration --method scvi --config config/atlas.yaml` | `processed/integrated_atlas.h5ad` |
| Thu | Annotation | `python -m src.atlas.annotate --config config/atlas.yaml` | `processed/integrated_annotated.h5ad` |
| Fri | Visualizations & proportions | `python -m src.atlas.viz --all --config config/atlas.yaml` | `processed/figures/*.png` |
| Sat | Export for cellxgene + Streamlit | `python -m src.atlas.export --cellxgene --config config/atlas.yaml` and `make app` | `processed/cellxgene_release/atlas.h5ad` |
| Sun | Final cleanup & docs | `make lint && make test && make all` | Updated README + repo push |

FAQ and troubleshooting
-----------------------

- Scanpy memory usage: Use backed mode for very large objects (`anndata.read_h5ad(..., backed='r')`) or subset genes/cells.
- Neighbors/UMAP determinism: We set `random_state` in neighbors and UMAP; minor variations can still occur.
- scVI with/without GPU: The pipeline will run on CPU by default. If GPU is available, install `pytorch-cuda` and scVI will detect it; otherwise pass `use_gpu=False`.
- Harmony vs scVI: Choose via `integration.method` in `config/atlas.yaml`. Harmony is faster and lighter; scVI can better model complex batch effects.
- Mitochondrial genes: Assumed to be named with the prefix `MT-` (human). Ensure gene names are in `var_names` and uppercased if needed.

Demo data
---------
The `make demo` target runs on synthetic data to keep the footprint small and ensure a quick end-to-end run. For real analyses, replace placeholder URLs in `config/atlas.yaml` with real datasets and re-run.

License
-------
This project is licensed under the MIT License (see `LICENSE`).

Contributing
------------
Pull requests are welcome. Please run `make lint` and `make test` before submitting changes.

TODOs for real-data runs
------------------------
- Replace placeholder dataset URLs in `config/atlas.yaml` with real URLs or local file paths
- Consider adding additional immune markers to `src/atlas/markers/immune_markers_human.tsv`
- Optionally enable GPU for scVI by installing `pytorch-cuda` (see `env.yml`)
