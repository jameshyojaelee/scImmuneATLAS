"""Cell type annotation functions for the Single-cell Immune Atlas."""

import logging
from pathlib import Path
from typing import Dict, List

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from .utils import load_config, setup_logging, timer


def load_immune_markers() -> pd.DataFrame:
    """Load immune cell type markers."""
    markers_path = Path(__file__).parent / "markers" / "immune_markers_human.tsv"
    
    if not markers_path.exists():
        raise FileNotFoundError(f"Markers file not found: {markers_path}")
    
    markers_df = pd.read_csv(markers_path, sep="\t")
    
    # Validate columns
    required_cols = ["cell_type", "gene", "direction"]
    if not all(col in markers_df.columns for col in required_cols):
        raise ValueError(f"Markers file must contain columns: {required_cols}")
    
    logging.info(f"Loaded {len(markers_df)} markers for {markers_df['cell_type'].nunique()} cell types")
    
    return markers_df


def compute_marker_scores(
    adata: ad.AnnData, 
    markers_df: pd.DataFrame, 
    min_pct: float = 0.1
) -> Dict[str, np.ndarray]:
    """Compute marker scores for each cell type."""
    scores = {}
    
    for cell_type in markers_df["cell_type"].unique():
        # Get markers for this cell type
        cell_markers = markers_df[markers_df["cell_type"] == cell_type]
        
        # Filter markers to those present in the dataset
        available_markers = cell_markers[cell_markers["gene"].isin(adata.var_names)]
        
        if len(available_markers) == 0:
            logging.warning(f"No markers found for {cell_type}")
            scores[cell_type] = np.zeros(adata.n_obs)
            continue
        
        # Compute mean expression of markers
        marker_genes = available_markers["gene"].tolist()
        marker_expr = adata[:, marker_genes].X
        
        if hasattr(marker_expr, "toarray"):
            marker_expr = marker_expr.toarray()
        
        # Simple scoring: mean expression
        cell_scores = np.mean(marker_expr, axis=1)
        scores[cell_type] = cell_scores
        
        logging.info(f"{cell_type}: {len(marker_genes)} markers, "
                    f"mean score: {np.mean(cell_scores):.3f}")
    
    return scores


def assign_cell_types(
    scores: Dict[str, np.ndarray], 
    min_score: float = 0.1
) -> pd.Series:
    """Assign cell types based on marker scores."""
    # Convert scores to DataFrame
    scores_df = pd.DataFrame(scores)
    
    # Find the cell type with highest score for each cell
    best_matches = scores_df.idxmax(axis=1)
    best_scores = scores_df.max(axis=1)
    
    # Filter out low-confidence assignments
    confident_mask = best_scores >= min_score
    
    # Assign "Unknown" to low-confidence cells
    cell_types = best_matches.copy()
    cell_types[~confident_mask] = "Unknown"
    
    # Log assignment statistics
    type_counts = cell_types.value_counts()
    logging.info("Cell type assignments:")
    for cell_type, count in type_counts.items():
        pct = 100 * count / len(cell_types)
        logging.info(f"  {cell_type}: {count} cells ({pct:.1f}%)")
    
    return cell_types


def score_and_annotate(adata: ad.AnnData, min_pct: float = 0.1) -> ad.AnnData:
    """Score markers and annotate cell types."""
    logging.info("Starting marker-based annotation")
    
    # Make a copy to avoid modifying original
    adata = adata.copy()
    
    # Load markers
    markers_df = load_immune_markers()
    
    # Compute marker scores
    scores = compute_marker_scores(adata, markers_df, min_pct)
    
    # Add scores to obsm
    for cell_type, score_array in scores.items():
        adata.obs[f"{cell_type}_score"] = score_array
    
    # Assign cell types
    cell_types = assign_cell_types(scores)
    adata.obs["cell_type"] = cell_types
    
    # Add broader compartment annotations
    compartment_map = {
        "CD8_T": "T_cell",
        "CD4_T": "T_cell", 
        "Treg": "T_cell",
        "NK": "NK_cell",
        "B_cell": "B_cell",
        "Plasma": "B_cell",
        "Mono": "Myeloid",
        "DC_cDC1": "Myeloid",
        "DC_cDC2": "Myeloid",
        "Unknown": "Unknown"
    }
    
    adata.obs["compartment"] = adata.obs["cell_type"].map(compartment_map)
    
    logging.info("Annotation completed")
    
    return adata


def run_annotation(config: dict) -> None:
    """Run cell type annotation."""
    # Load integrated data
    input_path = "processed/integrated_atlas.h5ad"
    logging.info(f"Loading integrated data from {input_path}")
    
    adata = ad.read_h5ad(input_path)
    
    # Run annotation
    adata_annotated = score_and_annotate(adata, config["annotation"]["min_pct"])
    
    # Save annotated data
    output_path = "processed/integrated_annotated.h5ad"
    adata_annotated.write(output_path)
    
    logging.info(f"Saved annotated atlas to {output_path}")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Run cell type annotation")
    parser.add_argument("--config", default="config/atlas.yaml", help="Config file path")
    
    args = parser.parse_args()
    
    setup_logging()
    config = load_config(args.config)
    
    with timer("Cell type annotation"):
        run_annotation(config)
