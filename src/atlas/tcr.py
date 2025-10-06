"""TCR repertoire analysis using Scirpy."""

from __future__ import annotations

import json
import logging
import tempfile
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .receptor.config import (
    dataset_has_receptor,
    dataset_receptor_config,
    figures_dir as config_figures_dir,
    global_receptor_config,
    metrics_dir as config_metrics_dir,
)
from .receptor.ingest import resolved_output_path as receptor_table_path
from .receptors_io import load_receptor_table
from .tcr_viz import generate_all_tcr_figures
from .utils import ensure_dir, load_config, set_seed, setup_logging, timer

logger = logging.getLogger(__name__)


def _import_scirpy():
    try:
        import scirpy as ir  # type: ignore
    except ImportError as exc:  # pragma: no cover - guarded import
        raise RuntimeError(
            "scirpy is required for TCR analysis. Please install scirpy within the analysis environment."
        ) from exc
    return ir


def _airr_from_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Convert harmonised receptor dataframe to AIRR-compliant schema for Scirpy."""
    airr_df = pd.DataFrame(
        {
            "cell_id": df["cell_id"],
            "locus": df.get("locus", df.get("chain", "")).astype(str),
            "junction": df.get("cdr3_nt", ""),
            "junction_aa": df.get("cdr3_aa", ""),
            "v_call": df.get("v_gene", ""),
            "d_call": df.get("d_gene", ""),
            "j_call": df.get("j_gene", ""),
            "productive": df.get("productive", pd.Series([pd.NA] * len(df))).astype(object),
            "duplicate_count": df.get("umis", pd.Series([pd.NA] * len(df))).astype(object),
            "consensus_count": df.get("reads", pd.Series([pd.NA] * len(df))).astype(object),
            "clonotype_id": df.get("clonotype_id", ""),
        }
    )
    return airr_df


def _merge_with_scirpy(adata: ad.AnnData, receptor_df: pd.DataFrame) -> ad.AnnData:
    ir = _import_scirpy()
    airr_df = _airr_from_dataframe(receptor_df)
    with tempfile.TemporaryDirectory() as tmpdir:
        airr_path = Path(tmpdir) / "receptor_airr.tsv"
        airr_df.to_csv(airr_path, sep="\t", index=False)
        ir.pp.merge_with_ir(adata, airr_path, how="left")
    return adata


def _compute_dataset_metrics(
    dataset_id: str,
    adata: ad.AnnData,
    receptor_df: pd.DataFrame,
    *,
    config: Dict,
) -> Tuple[Dict, Dict]:
    """Compute TCR metrics and return both metrics dict and diversity dict for visualization.

    Returns:
        Tuple of (metrics_dict, diversity_dict)
    """
    ir = _import_scirpy()

    # Define clonotypes & compute expansion/diversity
    ir.tl.define_clonotypes(adata)
    try:
        expansion = ir.tl.clonal_expansion(adata, inplace=False)
    except TypeError:  # pragma: no cover - API shim
        expansion = ir.tl.clonal_expansion(adata, target_col="clonet", inplace=False)

    expansion_counts = (
        expansion if isinstance(expansion, pd.Series) else pd.Series(expansion)
    ).fillna(0).astype(int)
    expansion_summary = expansion_counts.to_dict()

    diversity_df = ir.tl.alpha_diversity(adata, groupby=None, target_col="clonotype")
    if isinstance(diversity_df, pd.DataFrame):
        diversity = diversity_df.iloc[0].to_dict()
    else:  # pragma: no cover - Series fallback
        diversity = diversity_df.to_dict()

    # V/J usage using receptor dataframe for transparency
    v_usage = (
        receptor_df["v_gene"].replace("", pd.NA).dropna().value_counts().head(20).to_dict()
    )
    j_usage = (
        receptor_df["j_gene"].replace("", pd.NA).dropna().value_counts().head(20).to_dict()
    )

    clonotype_counts = (
        adata.obs.get("clonotype")
        .replace(["", None], pd.NA)
        .dropna()
        .value_counts()
        .to_dict()
    )

    metrics = {
        "dataset_id": dataset_id,
        "n_cells": int(adata.n_obs),
        "n_clonotypes": int(len(clonotype_counts)),
        "clonal_expansion": expansion_summary,
        "diversity": {
            key: float(value) for key, value in diversity.items() if pd.notna(value)
        },
        "v_gene_usage": v_usage,
        "j_gene_usage": j_usage,
        "top_clonotypes": sorted(
            (
                {"clonotype": clonotype, "cells": int(count)}
                for clonotype, count in clonotype_counts.items()
            ),
            key=lambda item: item["cells"],
            reverse=True,
        )[:10],
    }

    # Return diversity dict for visualization aggregation
    diversity_for_viz = {
        "dataset_id": dataset_id,
        **diversity,
    }

    return metrics, diversity_for_viz


def _load_dataset_artifacts(dataset_cfg: Dict, config: Dict) -> Tuple[ad.AnnData, pd.DataFrame]:
    dataset_id = dataset_cfg["id"]
    adata_path = Path("data/interim") / f"{dataset_id}.doublet_filtered.h5ad"
    if not adata_path.exists():
        raise FileNotFoundError(
            f"Doublet-filtered AnnData not found for dataset '{dataset_id}': {adata_path}"
        )
    adata = ad.read_h5ad(adata_path)

    receptor_table = receptor_table_path(dataset_cfg, config)
    if not receptor_table.exists():
        raise FileNotFoundError(
            f"Receptor table missing for dataset '{dataset_id}': {receptor_table}"
        )
    receptor_df = pd.read_parquet(receptor_table)
    return adata, receptor_df


def _write_json(path: Path, payload: Dict) -> None:
    ensure_dir(path.parent)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True))


def _calculate_jaccard_index(set1: set, set2: set) -> float:
    """Calculate Jaccard similarity between two sets."""
    if not set1 and not set2:
        return 0.0
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    return intersection / union if union > 0 else 0.0


def _calculate_morisita_horn(counts1: Counter, counts2: Counter) -> float:
    """Calculate Morisita-Horn similarity index between two abundance distributions."""
    all_species = set(counts1.keys()) | set(counts2.keys())
    if not all_species:
        return 0.0

    n1 = sum(counts1.values())
    n2 = sum(counts2.values())

    if n1 == 0 or n2 == 0:
        return 0.0

    sum_product = sum(counts1.get(sp, 0) * counts2.get(sp, 0) for sp in all_species)
    sum_sq1 = sum(v * v for v in counts1.values())
    sum_sq2 = sum(v * v for v in counts2.values())

    lambda1 = sum_sq1 / (n1 * n1)
    lambda2 = sum_sq2 / (n2 * n2)

    denominator = (lambda1 + lambda2) * n1 * n2
    if denominator == 0:
        return 0.0

    return (2 * sum_product) / denominator


def _compute_repertoire_overlap(
    datasets_receptor: Dict[str, pd.DataFrame]
) -> Dict:
    """Compute pairwise repertoire overlap metrics between datasets."""
    dataset_ids = list(datasets_receptor.keys())
    n_datasets = len(dataset_ids)

    # Extract clonotypes per dataset
    dataset_clonotypes = {}
    dataset_counts = {}

    for dataset_id, df in datasets_receptor.items():
        clonotypes = df["clonotype_id"].replace("", pd.NA).dropna()
        dataset_clonotypes[dataset_id] = set(clonotypes)
        dataset_counts[dataset_id] = Counter(clonotypes)

    # Compute pairwise metrics
    jaccard_matrix = np.zeros((n_datasets, n_datasets))
    morisita_matrix = np.zeros((n_datasets, n_datasets))

    for i, ds1 in enumerate(dataset_ids):
        for j, ds2 in enumerate(dataset_ids):
            if i == j:
                jaccard_matrix[i, j] = 1.0
                morisita_matrix[i, j] = 1.0
            else:
                jaccard_matrix[i, j] = _calculate_jaccard_index(
                    dataset_clonotypes[ds1], dataset_clonotypes[ds2]
                )
                morisita_matrix[i, j] = _calculate_morisita_horn(
                    dataset_counts[ds1], dataset_counts[ds2]
                )

    return {
        "dataset_ids": dataset_ids,
        "jaccard": jaccard_matrix.tolist(),
        "morisita_horn": morisita_matrix.tolist(),
    }


def _identify_public_clonotypes(
    datasets_receptor: Dict[str, pd.DataFrame],
    min_datasets: int = 2,
) -> Dict:
    """Identify clonotypes shared across multiple datasets (public clonotypes)."""
    # Collect all clonotypes with their dataset origins
    clonotype_datasets = {}
    clonotype_cell_counts = Counter()

    for dataset_id, df in datasets_receptor.items():
        clonotypes = df["clonotype_id"].replace("", pd.NA).dropna()
        for clonotype in clonotypes.unique():
            if clonotype not in clonotype_datasets:
                clonotype_datasets[clonotype] = set()
            clonotype_datasets[clonotype].add(dataset_id)
            clonotype_cell_counts[clonotype] += (clonotypes == clonotype).sum()

    # Filter for public clonotypes
    public_clonotypes = {
        clono: sorted(datasets)
        for clono, datasets in clonotype_datasets.items()
        if len(datasets) >= min_datasets
    }

    # Sort by number of datasets and cell count
    public_sorted = sorted(
        public_clonotypes.items(),
        key=lambda x: (len(x[1]), clonotype_cell_counts[x[0]]),
        reverse=True,
    )

    return {
        "n_public_clonotypes": len(public_clonotypes),
        "min_datasets_threshold": min_datasets,
        "top_public_clonotypes": [
            {
                "clonotype_id": clono,
                "n_datasets": len(datasets),
                "datasets": datasets,
                "total_cells": int(clonotype_cell_counts[clono]),
            }
            for clono, datasets in public_sorted[:50]
        ],
        "sharing_distribution": dict(
            Counter(len(datasets) for datasets in public_clonotypes.values())
        ),
    }


def _analyze_cdr3_properties(receptor_df: pd.DataFrame) -> Dict:
    """Analyze physicochemical properties of CDR3 sequences."""
    # Amino acid properties
    aa_properties = {
        'charge': {
            'R': 1, 'K': 1, 'D': -1, 'E': -1, 'H': 0.5,
            'A': 0, 'C': 0, 'F': 0, 'G': 0, 'I': 0, 'L': 0,
            'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'S': 0, 'T': 0,
            'V': 0, 'W': 0, 'Y': 0,
        },
        'hydrophobicity': {  # Kyte-Doolittle scale
            'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5,
            'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'S': -0.8,
            'W': -0.9, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5,
            'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5,
        },
    }

    cdr3_aa = receptor_df["cdr3_aa"].replace("", pd.NA).dropna()

    if cdr3_aa.empty:
        return {
            "mean_length": 0,
            "median_length": 0,
            "mean_charge": 0,
            "mean_hydrophobicity": 0,
        }

    lengths = cdr3_aa.str.len()

    # Calculate charge
    charges = []
    hydrophobicities = []

    for seq in cdr3_aa:
        if isinstance(seq, str) and seq:
            charge = sum(aa_properties['charge'].get(aa, 0) for aa in seq)
            charges.append(charge / len(seq))  # normalized by length

            hydro = sum(aa_properties['hydrophobicity'].get(aa, 0) for aa in seq)
            hydrophobicities.append(hydro / len(seq))

    return {
        "mean_length": float(lengths.mean()),
        "median_length": float(lengths.median()),
        "std_length": float(lengths.std()),
        "mean_charge": float(np.mean(charges)) if charges else 0.0,
        "mean_hydrophobicity": float(np.mean(hydrophobicities)) if hydrophobicities else 0.0,
        "length_distribution": dict(lengths.value_counts().head(20).to_dict()),
    }


def _generate_clonotype_network(
    adata: ad.AnnData,
    figures_dir: Path,
    max_clonotypes: int = 100,
) -> None:
    """Generate clonotype network visualization using scirpy."""
    ir = _import_scirpy()

    try:
        # Define clonotype networks based on sequence similarity
        # Only use most abundant clonotypes to avoid overcrowding
        clonotype_counts = adata.obs["clonotype"].value_counts()
        top_clonotypes = clonotype_counts.head(max_clonotypes).index.tolist()

        # Filter AnnData to top clonotypes
        mask = adata.obs["clonotype"].isin(top_clonotypes)
        if mask.sum() < 10:  # Need minimum cells for meaningful network
            logger.warning("Too few cells with top clonotypes for network analysis")
            return

        adata_subset = adata[mask].copy()

        # Compute clonotype network
        ir.tl.define_clonotype_clusters(
            adata_subset,
            sequence="junction_aa",
            metric="alignment",
            receptor_arms="all",
        )

        # Generate network plot
        fig, ax = plt.subplots(figsize=(12, 10))
        try:
            ir.pl.clonotype_network(
                adata_subset,
                color="clonotype_cluster",
                base_size=5,
                show=False,
                ax=ax,
            )
        except TypeError:  # API compatibility
            ir.pl.clonotype_network(
                adata_subset,
                color="clonotype_cluster",
                show=False,
                ax=ax,
            )

        ax.set_title(f"Clonotype Similarity Network (Top {max_clonotypes} Clonotypes)")
        plt.tight_layout()
        plt.savefig(figures_dir / "clonotype_network.png", dpi=300, bbox_inches="tight")
        plt.close()

        logger.info("Generated clonotype network visualization")

    except Exception as e:
        logger.warning("Failed to generate clonotype network: %s", e)
        # Create placeholder
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.axis("off")
        ax.text(
            0.5, 0.5,
            f"Clonotype network generation failed:\n{str(e)}",
            ha="center", va="center", wrap=True
        )
        plt.savefig(figures_dir / "clonotype_network.png", dpi=300, bbox_inches="tight")
        plt.close()


def _generate_tcr_visualizations(
    datasets_receptor: Dict[str, pd.DataFrame],
    overlap_metrics: Dict,
    public_clonotypes: Dict,
    figures_dir: Path,
    adata_dict: Optional[Dict[str, ad.AnnData]] = None,
) -> None:
    """Generate TCR-specific visualization figures."""
    ensure_dir(figures_dir)

    # 1. Repertoire overlap heatmap
    dataset_ids = overlap_metrics["dataset_ids"]
    jaccard_matrix = np.array(overlap_metrics["jaccard"])

    fig, ax = plt.subplots(figsize=(8, 7))
    sns.heatmap(
        jaccard_matrix,
        annot=True,
        fmt=".3f",
        cmap="viridis",
        xticklabels=dataset_ids,
        yticklabels=dataset_ids,
        square=True,
        cbar_kws={"label": "Jaccard Index"},
        ax=ax,
    )
    ax.set_title("Repertoire Overlap (Jaccard Index)")
    plt.tight_layout()
    plt.savefig(figures_dir / "repertoire_overlap_jaccard.png", dpi=300, bbox_inches="tight")
    plt.close()

    # 2. Morisita-Horn similarity heatmap
    morisita_matrix = np.array(overlap_metrics["morisita_horn"])

    fig, ax = plt.subplots(figsize=(8, 7))
    sns.heatmap(
        morisita_matrix,
        annot=True,
        fmt=".3f",
        cmap="plasma",
        xticklabels=dataset_ids,
        yticklabels=dataset_ids,
        square=True,
        cbar_kws={"label": "Morisita-Horn Index"},
        ax=ax,
    )
    ax.set_title("Repertoire Similarity (Morisita-Horn Index)")
    plt.tight_layout()
    plt.savefig(figures_dir / "repertoire_similarity_morisita.png", dpi=300, bbox_inches="tight")
    plt.close()

    # 3. Public clonotype distribution
    if public_clonotypes["n_public_clonotypes"] > 0:
        sharing_dist = public_clonotypes["sharing_distribution"]

        fig, ax = plt.subplots(figsize=(8, 5))
        x = sorted(sharing_dist.keys())
        y = [sharing_dist[k] for k in x]
        ax.bar(x, y, color="steelblue", edgecolor="black")
        ax.set_xlabel("Number of Datasets Sharing Clonotype")
        ax.set_ylabel("Number of Clonotypes")
        ax.set_title("Public Clonotype Sharing Distribution")
        ax.set_xticks(x)
        plt.tight_layout()
        plt.savefig(figures_dir / "public_clonotype_distribution.png", dpi=300, bbox_inches="tight")
        plt.close()

        # 4. Top public clonotypes
        top_public = pd.DataFrame(public_clonotypes["top_public_clonotypes"][:20])
        if not top_public.empty:
            fig, ax = plt.subplots(figsize=(10, 6))
            sns.barplot(
                data=top_public,
                x="n_datasets",
                y="clonotype_id",
                hue="total_cells",
                palette="viridis",
                dodge=False,
                ax=ax,
            )
            ax.set_xlabel("Number of Datasets")
            ax.set_ylabel("Clonotype ID")
            ax.set_title("Top Public Clonotypes")
            ax.legend(title="Total Cells", loc="lower right")
            plt.tight_layout()
            plt.savefig(figures_dir / "top_public_clonotypes.png", dpi=300, bbox_inches="tight")
            plt.close()

    # 5. CDR3 length distribution across datasets
    fig, ax = plt.subplots(figsize=(10, 6))
    for dataset_id, df in datasets_receptor.items():
        cdr3_lengths = df["cdr3_aa"].replace("", pd.NA).dropna().str.len()
        if not cdr3_lengths.empty:
            ax.hist(cdr3_lengths, bins=30, alpha=0.5, label=dataset_id, density=True)

    ax.set_xlabel("CDR3 Length (AA)")
    ax.set_ylabel("Density")
    ax.set_title("CDR3 Length Distribution by Dataset")
    ax.legend()
    plt.tight_layout()
    plt.savefig(figures_dir / "cdr3_length_distribution.png", dpi=300, bbox_inches="tight")
    plt.close()

    # 6. Clonotype network (if AnnData objects provided)
    if adata_dict:
        # Use the largest dataset for network visualization
        largest_dataset = max(adata_dict.items(), key=lambda x: x[1].n_obs)
        logger.info("Generating clonotype network from dataset: %s", largest_dataset[0])
        _generate_clonotype_network(largest_dataset[1], figures_dir)

    logger.info("Generated TCR visualization figures in %s", figures_dir)


def run_tcr_analysis(config: Dict) -> Dict:
    """Execute scirpy-based TCR analysis across datasets."""

    receptor_cfg = global_receptor_config(config)
    if not receptor_cfg.get("enabled", False):
        logger.info("TCR analysis disabled in configuration; skipping")
        return {}

    set_seed(config.get("seed", 0))
    metrics_dir = config_metrics_dir(config)
    figures_dir = config_figures_dir(config)
    ensure_dir(metrics_dir)
    ensure_dir(figures_dir)

    datasets = [
        dataset_cfg
        for dataset_cfg in config.get("datasets", [])
        if dataset_has_receptor(dataset_cfg)
    ]

    if not datasets:
        logger.warning("TCR analysis requested but no datasets provide receptor data")
        return {}

    dataset_metrics: Dict[str, Dict] = {}
    aggregated_receptor: List[pd.DataFrame] = []
    datasets_receptor: Dict[str, pd.DataFrame] = {}
    adata_dict: Dict[str, ad.AnnData] = {}
    diversity_records: List[Dict] = []

    for dataset_cfg in datasets:
        dataset_id = dataset_cfg["id"]
        try:
            with timer(f"TCR analysis ({dataset_id})"):
                adata, receptor_df = _load_dataset_artifacts(dataset_cfg, config)
                aggregated_receptor.append(receptor_df.assign(dataset_id=dataset_id))
                datasets_receptor[dataset_id] = receptor_df
                _merge_with_scirpy(adata, receptor_df)
                adata_dict[dataset_id] = adata  # Store for network visualization
                metrics, diversity_dict = _compute_dataset_metrics(dataset_id, adata, receptor_df, config=config)

                # Add CDR3 property analysis per dataset
                cdr3_props = _analyze_cdr3_properties(receptor_df)
                metrics["cdr3_properties"] = cdr3_props

                # Store diversity metrics for visualization
                diversity_records.append(diversity_dict)

                dataset_metrics[dataset_id] = metrics
                _write_json(metrics_dir / f"{dataset_id}_tcr_metrics.json", metrics)
        except Exception as exc:  # pragma: no cover - runtime safety
            logger.exception("Failed TCR analysis for %s: %s", dataset_id, exc)
            raise

    aggregated_df = pd.concat(aggregated_receptor, ignore_index=True)

    # Compute repertoire overlap metrics
    with timer("Repertoire overlap analysis"):
        overlap_metrics = _compute_repertoire_overlap(datasets_receptor)
        _write_json(metrics_dir / "repertoire_overlap.json", overlap_metrics)
        logger.info("Computed repertoire overlap for %d datasets", len(datasets_receptor))

    # Identify public clonotypes
    with timer("Public clonotype detection"):
        min_datasets = receptor_cfg.get("min_public_datasets", 2)
        public_clonotypes = _identify_public_clonotypes(datasets_receptor, min_datasets=min_datasets)
        _write_json(metrics_dir / "public_clonotypes.json", public_clonotypes)
        logger.info(
            "Identified %d public clonotypes (shared across ≥%d datasets)",
            public_clonotypes["n_public_clonotypes"],
            min_datasets,
        )

    # Global CDR3 properties
    with timer("CDR3 property analysis"):
        global_cdr3_props = _analyze_cdr3_properties(aggregated_df)

    global_metrics = {
        "total_cells": int(sum(m["n_cells"] for m in dataset_metrics.values())),
        "total_clonotypes": int(
            aggregated_df["clonotype_id"].replace("", pd.NA).dropna().nunique()
        ),
        "v_gene_usage": (
            aggregated_df["v_gene"].replace("", pd.NA).dropna().value_counts().head(25).to_dict()
        ),
        "j_gene_usage": (
            aggregated_df["j_gene"].replace("", pd.NA).dropna().value_counts().head(25).to_dict()
        ),
        "cdr3_properties": global_cdr3_props,
        "n_public_clonotypes": public_clonotypes["n_public_clonotypes"],
        "repertoire_overlap_summary": {
            "mean_jaccard": float(np.mean([
                overlap_metrics["jaccard"][i][j]
                for i in range(len(overlap_metrics["dataset_ids"]))
                for j in range(i + 1, len(overlap_metrics["dataset_ids"]))
            ])) if len(overlap_metrics["dataset_ids"]) > 1 else 0.0,
            "mean_morisita_horn": float(np.mean([
                overlap_metrics["morisita_horn"][i][j]
                for i in range(len(overlap_metrics["dataset_ids"]))
                for j in range(i + 1, len(overlap_metrics["dataset_ids"]))
            ])) if len(overlap_metrics["dataset_ids"]) > 1 else 0.0,
        },
    }

    # Generate visualizations
    with timer("TCR visualization"):
        # Original overlap/network visualizations
        _generate_tcr_visualizations(
            datasets_receptor,
            overlap_metrics,
            public_clonotypes,
            figures_dir,
            adata_dict=adata_dict,
        )

        # New comprehensive visualization suite
        try:
            # Try to load integrated annotated atlas for better UMAP with cancer_type
            atlas_path = Path(config.get("outputs", {}).get("integrated_path", "processed/integrated_annotated.h5ad"))
            groupby_col = "cancer_type"

            if atlas_path.exists():
                logger.info("Loading integrated atlas from %s for enhanced visualizations", atlas_path)
                atlas_adata = ad.read_h5ad(atlas_path)

                # Prepare diversity dataframe with grouping variable
                diversity_df = pd.DataFrame(diversity_records)
                if groupby_col in atlas_adata.obs.columns:
                    # Map dataset_id to cancer_type
                    dataset_to_cancer = atlas_adata.obs.groupby("dataset_id")[groupby_col].first().to_dict()
                    diversity_df[groupby_col] = diversity_df["dataset_id"].map(dataset_to_cancer)

                generate_all_tcr_figures(
                    atlas_adata,
                    aggregated_df,
                    diversity_df,
                    figures_dir,
                    groupby=groupby_col,
                )
            else:
                # Use largest individual dataset if no integrated atlas
                logger.info("Integrated atlas not found; using largest dataset for visualizations")
                largest_dataset = max(adata_dict.items(), key=lambda x: x[1].n_obs)
                diversity_df = pd.DataFrame(diversity_records)

                generate_all_tcr_figures(
                    largest_dataset[1],
                    aggregated_df,
                    diversity_df,
                    figures_dir,
                    groupby="dataset_id",
                )
        except Exception as e:
            logger.warning("Failed to generate comprehensive visualization suite: %s", e)
            logger.debug("Visualization error details:", exc_info=True)

    summary = {
        "datasets": dataset_metrics,
        "global": global_metrics,
        "overlap": overlap_metrics,
        "public_clonotypes": public_clonotypes,
    }
    _write_json(metrics_dir / "tcr_summary.json", summary)
    logger.info("TCR analysis complete for %d datasets", len(dataset_metrics))
    return summary


def main() -> None:  # pragma: no cover - CLI helper
    import argparse

    parser = argparse.ArgumentParser(description="Run TCR repertoire analysis with Scirpy")
    parser.add_argument("--config", default="config/atlas.yaml", help="YAML configuration path")
    parser.add_argument("--log-level", default="INFO", help="Logging level (default: INFO)")
    args = parser.parse_args()

    setup_logging(args.log_level)
    config = load_config(args.config)
    run_tcr_analysis(config)


if __name__ == "__main__":  # pragma: no cover
    main()
