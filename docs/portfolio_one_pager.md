# scImmuneATLAS – Interview One-Pager

## Snapshot
- **Scope:** Harmonised tumour-infiltrating immune atlas spanning melanoma, NSCLC, breast, and RCC cohorts.
- **Tech Stack:** Snakemake workflow, Scanpy/ScVI, Scrublet, Harmony, reproducible Conda environment.
- **Deliverables:** Annotated `.h5ad`, QC/doublet diagnostics, integration metrics, benchmarking report, interpretation notebook, cellxgene package.

## Quantitative Highlights
- QC retention: see `processed/metrics/*_qc_summary.json` (cell retention ~ TBD once run).
- Doublet removal: Scrublet histograms + summaries under `processed/metrics/doublets/`.
- Integration quality: silhouette + LISI metrics in `processed/metrics/integration_metrics.json`.
- Annotation coverage: `annotation_summary.json` (counts per immune lineage).

## Talking Points
- Marker scoring honours directionality and confidence gating; Unknown cells routed to reference fallback.
- Integration keeps sparse counts for scVI, selects shared HVGs, and outputs diagnostics.
- CLI (`scimmuneatlas <stage>`) standardises execution; Snakemake enforces metrics outputs for regression testing.
- Portfolio assets: interpretation notebook, clinical metadata scaffold, video/blog outlines ready for polish.

## Next Steps Before Demo
1. Run `make all` (or targeted CLI commands) to regenerate outputs with real data.
2. Populate `data/external/clinical_template.csv` with cohort metadata and refresh annotation confusion matrices.
3. Record video walkthrough following `docs/video_walkthrough_script.md` and publish via GitHub Pages / blog outline.
