# scImmuneATLAS – Interview One-Pager

## Snapshot
- **Scope:** Harmonised tumour-infiltrating immune atlas spanning melanoma, NSCLC, breast, and RCC cohorts.
- **Pipeline Toggles:** Harmony ↔ scVI integration switch, shared-HVG pre-selection, receptor/TCR analytics module (scirpy + custom summarisation).
- **Tech Stack:** Snakemake + typed CLI, Scanpy, scVI-tools, Harmonypy, Scrublet, Pandera schemas, reproducible Conda environment.
- **Deliverables:** Annotated `.h5ad`, QC/doublet diagnostics, integration metrics, receptor analytics, benchmarking-ready report, interpretation notebook, cellxgene package.

## Quantitative Highlights
- QC retention: see `processed/metrics/*_qc_summary.json` (cell retention ~ TBD once run).
- Doublet removal: Scrublet histograms + summaries under `processed/metrics/doublets/`.
- Integration quality: silhouette + LISI metrics in `processed/metrics/integration_metrics.json`.
- Annotation coverage: `annotation_summary.json` (counts per immune lineage).

## Talking Points
- Marker scoring honours directionality, applies dataset-level confidence gating, and can fall back to reference labels for ambiguous clusters.
- Integration keeps sparse counts for scVI, exposes `integration.method` switch (Harmony vs scVI), and emits regression-friendly diagnostics.
- Receptor/TCR module harmonises contig tables, computes diversity/overlap/public clonotypes, and writes report fragments and figures.
- CLI (`scimmuneatlas <stage>`) standardises execution; Snakemake + regression tests guarantee required metrics/figures before marking runs successful.
- Portfolio assets: interpretation notebook, clinical metadata scaffold, video/blog outlines, resume blurbs, and quick-start thumbnails.

## Next Steps Before Demo
1. Run `make all` (or targeted CLI commands) to regenerate outputs with real data.
2. Populate `data/external/clinical_template.csv` with cohort metadata and refresh annotation confusion matrices.
3. Record video walkthrough following `docs/video_walkthrough_script.md` and publish via GitHub Pages / blog outline.
