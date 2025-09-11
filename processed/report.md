# Single-cell Immune Atlas Report

## Configuration
- Project: single_cell_immune_atlas
- Organism: human
- Integration method: harmony
- Seed: 1337

## Datasets
- **MELANOMA_CENSUS**: melanoma (mixed)
- **NSCLC_CENSUS**: NSCLC (mixed)
- **BREAST_CENSUS**: Breast (mixed)
- **NSCLC_PEMBROLIZUMAB_CENSUS**: NSCLC (mixed)
- **NSCLC_NIVO_IPI_CENSUS**: NSCLC (mixed)
- **RCC_COMBINATION_CENSUS**: RCC (mixed)

## QC Parameters
- Min genes per cell: 200
- Max genes per cell: 6000
- Max mitochondrial %: 15

## Integration Parameters
- Method: harmony
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
