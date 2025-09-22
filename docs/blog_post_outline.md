# Blog Post Outline: Building a Reproducible Tumour Immune Atlas

## Abstract
- Two-sentence summary linking clinical motivation and technical solution.

## 1. Motivation
- Fragmented tumour immune datasets and cross-cohort comparability challenges.
- Requirements from translational collaborators (QC, harmonisation, annotation, explainability).

## 2. Pipeline Architecture
- Snakemake DAG + config-driven design.
- Modular stages: QC, doublets, integration (Harmony/scVI), annotation, benchmarking.
- CLI + CI enhancements for production readiness.

## 3. Scientific Highlights
- QC plots showcasing improved retention vs. historical pipelines.
- Integration metrics (silhouette/LISI) and benchmarking (ARI/NMI) vs. public references.
- Annotation upgrades: marker directionality, confidence gating, reference fallback.

## 4. Demo & Visualisation
- Interpretation notebook walkthrough (UMAPs, differential expression ideas).
- Cellxgene export + Streamlit app (link instructions).

## 5. Future Directions / Call for Collaboration
- Add spatial multi-omics, TCR/BCR integration, therapy response prediction.
- Invite readers to explore repo, share datasets, or discuss roles.
