"""Cell type annotation functions for the Single-cell Immune Atlas."""

import logging
from pathlib import Path
from typing import Dict, List

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

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
    min_pct: float = 0.1,
) -> Dict[str, np.ndarray]:
    """Compute marker scores for each cell type.

    Robust to var names being Ensembl IDs by matching markers against
    `var['feature_name']` when available; otherwise uses `var_names`.
    """
    scores: Dict[str, np.ndarray] = {}

    # Select the gene symbol series to match on
    if "feature_name" in adata.var.columns:
        symbol_series = adata.var["feature_name"].astype(str)
    elif getattr(adata, "raw", None) is not None and adata.raw is not None and "feature_name" in adata.raw.var.columns:
        # Map feature_name from raw (pre-HVG) to current var index
        raw_var = adata.raw.var
        # Build a mapping from raw var_names to feature_name
        raw_map = raw_var["feature_name"].astype(str)
        # Align to current var_names; missing keys fall back to current var_names
        mapped = []
        for v in adata.var_names.astype(str):
            mapped.append(raw_map.get(v, v))
        symbol_series = pd.Series(mapped, index=adata.var_names, dtype=str)
    else:
        symbol_series = pd.Series(adata.var_names, index=adata.var_names, dtype=str)

    for cell_type in markers_df["cell_type"].unique():
        cell_markers = markers_df[markers_df["cell_type"] == cell_type]
        marker_symbols = cell_markers["gene"].astype(str).tolist()

        # Find indices of marker symbols present in the dataset
        present_mask = symbol_series.isin(marker_symbols)
        marker_idx = np.where(present_mask.values)[0]

        if marker_idx.size == 0:
            logging.warning(f"No markers found for {cell_type}")
            scores[cell_type] = np.zeros(adata.n_obs)
            continue

        subset_X = adata[:, marker_idx].X
        if hasattr(subset_X, "toarray"):
            subset_X = subset_X.toarray()

        cell_scores = np.mean(subset_X, axis=1)
        scores[cell_type] = cell_scores

        logging.info(
            f"{cell_type}: {len(marker_idx)} markers, mean score: {np.mean(cell_scores):.3f}"
        )

    return scores


def _augment_var_with_feature_names_from_interim(adata: ad.AnnData) -> None:
    """If integrated AnnData lacks feature symbols, try to populate them
    by scanning interim H5ADs for a mapping of feature_id/soma_joinid → feature_name.
    """
    if "feature_name" in adata.var.columns:
        return
    interim_dir = Path("data/interim")
    if not interim_dir.exists():
        return
    mapping = {}
    for f in sorted(interim_dir.glob("*.h5ad")):
        try:
            ds = ad.read_h5ad(f)
            var = ds.var
            name_col = "feature_name" if "feature_name" in var.columns else None
            if name_col:
                # map by feature_id
                if "feature_id" in var.columns:
                    keys = var["feature_id"].astype(str)
                    names = var[name_col].astype(str)
                    mapping.update(dict(zip(keys, names)))
                # map by soma_joinid
                if "soma_joinid" in var.columns:
                    keys = var["soma_joinid"].astype(str)
                    names = var[name_col].astype(str)
                    mapping.update(dict(zip(keys, names)))
        except Exception:
            continue
    if mapping:
        adata.var["feature_name"] = [mapping.get(str(v), str(v)) for v in adata.var_names]


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

    # Try to augment feature_name if missing using interim mappings
    if "feature_name" not in adata.var.columns:
        _augment_var_with_feature_names_from_interim(adata)

    # If Census gene symbols are provided, keep a symbol view for matching
    if "feature_name" in adata.var.columns:
        adata.var["symbol"] = adata.var["feature_name"].astype(str)
    else:
        adata.var["symbol"] = np.array(adata.var_names, dtype=str)
    
    # Load markers
    markers_df = load_immune_markers()
    
    # Compute marker scores
    scores = compute_marker_scores(adata, markers_df, min_pct)
    
    # Add scores to obs
    for cell_type, score_array in scores.items():
        adata.obs[f"{cell_type}_score"] = score_array
    
    # Assign predicted cell types
    pred_cell_types = assign_cell_types(scores).astype(str)
    adata.obs["cell_type_pred"] = pd.Categorical(pred_cell_types.fillna("Unknown"))
    
    # Preserve existing labels if present; otherwise use predictions
    if "cell_type" not in adata.obs.columns or adata.obs["cell_type"].isna().all():
        adata.obs["cell_type"] = adata.obs["cell_type_pred"]
    
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
    
    adata.obs["compartment"] = adata.obs["cell_type"].map(compartment_map).astype(str).fillna("Unknown")
    adata.obs["compartment"] = pd.Categorical(adata.obs["compartment"])
    
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
