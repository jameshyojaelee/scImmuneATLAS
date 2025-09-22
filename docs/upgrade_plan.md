# scImmuneATLAS Upgrade Plan

## Configuration Extensions
- `qc.min_counts` / `qc.max_counts`: optional UMI thresholds evaluated alongside existing `min_genes` / `max_genes`.
- `qc.min_cells_per_gene`: configurable gene filter (default 3).
- `qc.qc_plots_dir`: output directory for per-dataset QC artifacts.
- `doublets.mask_high_mt_genes`: boolean; when true, temporarily removes top mitochondrial/ribosomal genes before Scrublet.
- `doublets.scrublet_kwargs`: dictionary forwarded to `scrub.scrub_doublets` (e.g., `min_counts`, `min_cells`, `min_gene_variability_pctl`, `n_prin_comps`).
- `doublets.save_plots`: toggle for writing Scrublet histograms / score distributions.
- `integration.use_gpu`: `"auto" | true | false` to control scVI device placement.
- `integration.n_hvg`: shared number of HVGs for Harmony and scVI (default 3000).
- `integration.shared_hvg`: boolean to reuse HVG list between methods.
- `integration.metrics`: list of integration metrics to compute (e.g., `kbET`, `lisi`, `silhouette`).
- `annotation.score_method`: `"scanpy_score_genes" | "aucell"` selection.
- `annotation.min_confidence_pct`: percentile of per-cell max score required for assignment.
- `annotation.reference_mapping`: optional block (method + reference paths) for fallback label transfer.
- `outputs.metrics_dir`: new directory (`processed/metrics/`) for JSON/TSV diagnostics.
- `report.embed_figures`: toggle controlling whether plots are embedded in `processed/report.md`.
- `benchmarking.enabled`: global flag for benchmarking routines.
- `benchmarking.reference_atlas`: path or identifier for public atlas used in comparisons.
- `benchmarking.metrics`: list of metrics (ARI, NMI, LISI, etc.).

## Artifact Inventory
- `data/interim/{dataset}_qc.tsv`: per-cell QC table before/after filtering.
- `processed/figures/qc/{dataset}_umi_violin.png`: UMI/genes violin plots per dataset.
- `processed/figures/qc/{dataset}_mt_violin.png`: mitochondrial percentage violin plots.
- `processed/metrics/qc_summary.tsv`: consolidated QC statistics.
- `processed/metrics/doublets/{dataset}_summary.json`: Scrublet summary (scores, threshold, counts).
- `processed/figures/doublets/{dataset}_scrublet_hist.png`: histogram of simulated vs observed scores.
- `processed/metrics/integration_metrics.json`: integration diagnostics.
- `processed/metrics/annotation_confusion.tsv`: confusion matrix vs. provided labels/reference.
- `processed/metrics/benchmarking.json`: results comparing against reference atlas.
- `docs/gallery.ipynb`: notebook showcasing benchmarking and biological interpretation.

## Workflow Touchpoints
1. QC module to retain `n_genes_by_counts`, compute new count-based filters, and export summaries/plots.
2. Doublet module to reuse a single Scrublet instance per dataset and persist thresholds/diagnostics.
3. Integration module to maintain sparse counts, shared HVG list, updated scVI setup, GPU autodetect, and metrics computation hooks.
4. Annotation module to respect marker direction, implement background-controlled scoring, confidence gating, and optional reference transfer.
5. Benchmarking utilities integrating scIB/LISI metrics and reference atlas comparisons; callable from Snakemake.
6. Report generator to consume new metric artifacts and embed figures.
7. Snakemake workflow to add checkpoints/guards and new outputs for QC/metrics/benchmarking.
8. Packaging updates (CLI entry points, CI, documentation) to be tackled after core scientific improvements.

## Immediate Next Steps
- Update `config/atlas.yaml` with new defaults (non-breaking fallbacks where possible).
- Implement QC module changes and ensure tests cover new behavior.
- Introduce artifact directories (`processed/metrics`, `processed/figures/qc`, `processed/figures/doublets`).
