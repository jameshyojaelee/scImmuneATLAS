# Single-cell Immune Atlas Report

## Configuration
- Project: single_cell_immune_atlas
- Organism: human
- Integration method: scvi
- Seed: 1337

## Datasets
- **GSE98638**: melanoma (10x)
- **GSE131907**: NSCLC (10x)

## QC Parameters
- Min genes per cell: 200
- Max genes per cell: 6000
- Max mitochondrial %: 15

## Integration Parameters
- Method: scvi
- Latent dimensions: 30
- Batch key: dataset_id

## Outputs
- Integrated atlas: `processed/integrated_atlas.h5ad`
- Annotated atlas: `processed/integrated_annotated.h5ad`
- Cellxgene export: `processed/cellxgene_release/atlas.h5ad`
- Figures: `processed/figures/`

## Generated Figures
- UMAP by cell type
- UMAP by dataset
- UMAP by cancer type
- Cell type proportions by cancer type
