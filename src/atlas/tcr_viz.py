"""Visualization suite for TCR repertoire analysis.

This module provides specialized visualization functions for TCR/BCR analysis,
generating publication-quality figures with consistent styling.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .utils import ensure_dir

logger = logging.getLogger(__name__)

# Consistent styling
sns.set_style("whitegrid")
FIGURE_DPI = 300


def plot_clonotype_frequency(
    adata: ad.AnnData,
    clonotype_col: str = "clonotype",
    top_n: int = 20,
    outpath: Optional[Path] = None,
    dataset_col: Optional[str] = None,
    figsize: Tuple[int, int] = (10, 8),
) -> None:
    """Generate barplot of top-N expanded clonotypes.

    Args:
        adata: AnnData object with clonotype annotations
        clonotype_col: Column name containing clonotype IDs
        top_n: Number of top clonotypes to display
        outpath: Output file path for PNG
        dataset_col: Optional column for dataset grouping
        figsize: Figure size tuple
    """
    if clonotype_col not in adata.obs.columns:
        logger.warning("Clonotype column '%s' not found; skipping frequency plot", clonotype_col)
        return

    # Count clonotype frequencies
    clonotype_counts = adata.obs[clonotype_col].value_counts().head(top_n)

    if clonotype_counts.empty:
        logger.warning("No clonotypes found for frequency plot")
        return

    # Prepare data
    plot_data = pd.DataFrame({
        "clonotype": clonotype_counts.index,
        "cells": clonotype_counts.values,
    })

    # Add dataset information if available
    if dataset_col and dataset_col in adata.obs.columns:
        dataset_info = []
        for clono in plot_data["clonotype"]:
            mask = adata.obs[clonotype_col] == clono
            datasets = adata.obs.loc[mask, dataset_col].value_counts()
            dominant_dataset = datasets.index[0] if len(datasets) > 0 else "Unknown"
            dataset_info.append(dominant_dataset)
        plot_data["dataset"] = dataset_info

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Reverse order for horizontal bar plot (largest at top)
    plot_data = plot_data.iloc[::-1]

    if dataset_col and "dataset" in plot_data.columns:
        sns.barplot(
            data=plot_data,
            y="clonotype",
            x="cells",
            hue="dataset",
            dodge=False,
            palette="tab10",
            ax=ax,
        )
        ax.legend(title="Dataset", bbox_to_anchor=(1.05, 1), loc="upper left")
    else:
        sns.barplot(
            data=plot_data,
            y="clonotype",
            x="cells",
            color="steelblue",
            ax=ax,
        )

    ax.set_xlabel("Number of Cells", fontsize=12)
    ax.set_ylabel("Clonotype ID", fontsize=12)
    ax.set_title(f"Top {top_n} Expanded Clonotypes", fontsize=14, fontweight="bold")

    plt.tight_layout()

    if outpath:
        ensure_dir(outpath.parent)
        plt.savefig(outpath, dpi=FIGURE_DPI, bbox_inches="tight")
        logger.info("Saved clonotype frequency plot to %s", outpath)

    plt.close(fig)


def plot_repertoire_diversity(
    diversity_df: pd.DataFrame,
    groupby: str = "cancer_type",
    outpath: Optional[Path] = None,
    figsize: Tuple[int, int] = (12, 6),
) -> None:
    """Generate boxplots of repertoire diversity metrics per group.

    Args:
        diversity_df: DataFrame with columns [group, metric_name, metric_value]
        groupby: Column name for grouping (e.g., cancer_type, dataset_id)
        outpath: Output file path for PNG
        figsize: Figure size tuple
    """
    if diversity_df.empty:
        logger.warning("Empty diversity dataframe; skipping diversity plot")
        return

    # Identify available diversity metrics
    metric_cols = [col for col in diversity_df.columns if col not in [groupby, "dataset_id"]]

    if not metric_cols:
        logger.warning("No diversity metrics found in dataframe")
        return

    # Create subplots for each metric
    n_metrics = len(metric_cols)
    fig, axes = plt.subplots(1, n_metrics, figsize=(figsize[0], figsize[1]), squeeze=False)
    axes = axes.flatten()

    for idx, metric in enumerate(metric_cols):
        ax = axes[idx]

        # Filter valid values
        plot_data = diversity_df[[groupby, metric]].dropna()

        if plot_data.empty:
            ax.text(0.5, 0.5, f"No data for {metric}", ha="center", va="center")
            ax.set_xlabel(groupby)
            ax.set_ylabel(metric)
            continue

        # Create boxplot
        sns.boxplot(
            data=plot_data,
            x=groupby,
            y=metric,
            palette="Set2",
            ax=ax,
        )

        # Add individual points
        sns.stripplot(
            data=plot_data,
            x=groupby,
            y=metric,
            color="black",
            alpha=0.3,
            size=3,
            ax=ax,
        )

        ax.set_xlabel(groupby.replace("_", " ").title(), fontsize=11)
        ax.set_ylabel(metric.replace("_", " ").title(), fontsize=11)
        ax.set_title(metric.replace("_", " ").title(), fontsize=12, fontweight="bold")

        # Rotate x-axis labels if needed
        if len(plot_data[groupby].unique()) > 3:
            ax.tick_params(axis="x", rotation=45)

    plt.suptitle(f"Repertoire Diversity by {groupby.replace('_', ' ').title()}",
                 fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()

    if outpath:
        ensure_dir(outpath.parent)
        plt.savefig(outpath, dpi=FIGURE_DPI, bbox_inches="tight")
        logger.info("Saved repertoire diversity plot to %s", outpath)

    plt.close(fig)


def plot_umap_clonotype_size(
    adata: ad.AnnData,
    clonotype_col: str = "clonotype",
    umap_key: str = "X_umap",
    outpath: Optional[Path] = None,
    size_bins: List[int] = [1, 2, 5, 10, 50],
    figsize: Tuple[int, int] = (10, 8),
    max_points: int = 50000,
) -> None:
    """Generate UMAP colored by clonotype expansion size.

    Args:
        adata: AnnData object with UMAP embeddings
        clonotype_col: Column name containing clonotype IDs
        umap_key: Key in adata.obsm for UMAP coordinates
        outpath: Output file path for PNG
        size_bins: Bin edges for clonotype size categories
        figsize: Figure size tuple
        max_points: Maximum points to plot (downsample if needed)
    """
    if umap_key not in adata.obsm:
        logger.warning("UMAP coordinates '%s' not found; skipping UMAP plot", umap_key)
        return

    if clonotype_col not in adata.obs.columns:
        logger.warning("Clonotype column '%s' not found; skipping UMAP plot", clonotype_col)
        return

    # Get UMAP coordinates
    umap_coords = adata.obsm[umap_key]

    # Calculate clonotype sizes
    clonotype_sizes = adata.obs[clonotype_col].map(
        adata.obs[clonotype_col].value_counts()
    ).fillna(1).astype(int)

    # Create size bins
    size_labels = [f"Singleton"]
    size_labels.extend([f"{size_bins[i]}-{size_bins[i+1]-1}"
                        for i in range(len(size_bins)-1)])
    size_labels.append(f"≥{size_bins[-1]}")

    size_categories = pd.cut(
        clonotype_sizes,
        bins=[0] + size_bins + [float("inf")],
        labels=size_labels,
        right=False,
    )

    # Prepare plot dataframe
    plot_df = pd.DataFrame({
        "UMAP1": umap_coords[:, 0],
        "UMAP2": umap_coords[:, 1],
        "size_category": size_categories,
        "size_value": clonotype_sizes,
    })

    # Downsample if needed
    if len(plot_df) > max_points:
        logger.info("Downsampling from %d to %d points for UMAP plot", len(plot_df), max_points)
        plot_df = plot_df.sample(n=max_points, random_state=42)

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Plot 1: Categorical size bins
    for category in size_labels:
        mask = plot_df["size_category"] == category
        subset = plot_df[mask]
        ax1.scatter(
            subset["UMAP1"],
            subset["UMAP2"],
            s=10,
            alpha=0.6,
            label=category,
            rasterized=True,
        )

    ax1.set_xlabel("UMAP 1", fontsize=12)
    ax1.set_ylabel("UMAP 2", fontsize=12)
    ax1.set_title("Clonal Expansion (Categorical)", fontsize=12, fontweight="bold")
    ax1.legend(title="Clone Size", bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=9)
    ax1.set_aspect("equal")

    # Plot 2: Continuous size with colormap
    scatter = ax2.scatter(
        plot_df["UMAP1"],
        plot_df["UMAP2"],
        c=plot_df["size_value"],
        s=10,
        alpha=0.6,
        cmap="viridis",
        norm=plt.Normalize(vmin=1, vmax=min(plot_df["size_value"].quantile(0.95), 100)),
        rasterized=True,
    )

    ax2.set_xlabel("UMAP 1", fontsize=12)
    ax2.set_ylabel("UMAP 2", fontsize=12)
    ax2.set_title("Clonal Expansion (Continuous)", fontsize=12, fontweight="bold")
    ax2.set_aspect("equal")

    cbar = plt.colorbar(scatter, ax=ax2)
    cbar.set_label("Cells per Clonotype", fontsize=10)

    plt.suptitle("UMAP by Clonotype Size", fontsize=14, fontweight="bold", y=0.98)
    plt.tight_layout()

    if outpath:
        ensure_dir(outpath.parent)
        plt.savefig(outpath, dpi=FIGURE_DPI, bbox_inches="tight")
        logger.info("Saved UMAP clonotype size plot to %s", outpath)

    plt.close(fig)


def plot_spectratype(
    receptor_df: pd.DataFrame,
    cdr3_col: str = "cdr3_nt",
    chain_col: str = "chain",
    dataset_col: Optional[str] = "dataset_id",
    outpath: Optional[Path] = None,
    figsize: Tuple[int, int] = (12, 6),
    separate_chains: bool = True,
) -> None:
    """Generate CDR3 length distribution (spectratype) plots.

    Args:
        receptor_df: DataFrame with receptor data
        cdr3_col: Column name for CDR3 sequences (nucleotide)
        chain_col: Column name for chain type
        dataset_col: Optional column for dataset grouping
        outpath: Output file path for PNG
        figsize: Figure size tuple
        separate_chains: Whether to separate by chain type
    """
    if cdr3_col not in receptor_df.columns:
        logger.warning("CDR3 column '%s' not found; skipping spectratype plot", cdr3_col)
        return

    # Calculate CDR3 lengths
    cdr3_lengths = receptor_df[cdr3_col].str.len()
    plot_data = receptor_df.copy()
    plot_data["cdr3_length"] = cdr3_lengths

    # Remove invalid lengths
    plot_data = plot_data[plot_data["cdr3_length"].notna() & (plot_data["cdr3_length"] > 0)]

    if plot_data.empty:
        logger.warning("No valid CDR3 lengths found for spectratype plot")
        return

    # Determine layout
    if separate_chains and chain_col in plot_data.columns:
        chains = plot_data[chain_col].dropna().unique()
        chains = sorted([c for c in chains if c])

        if len(chains) == 0:
            separate_chains = False
        else:
            n_plots = len(chains)
            fig, axes = plt.subplots(1, n_plots, figsize=figsize, squeeze=False)
            axes = axes.flatten()

            for idx, chain in enumerate(chains):
                ax = axes[idx]
                chain_data = plot_data[plot_data[chain_col] == chain]

                if dataset_col and dataset_col in chain_data.columns:
                    # Plot by dataset
                    for dataset in chain_data[dataset_col].unique():
                        dataset_lengths = chain_data[chain_data[dataset_col] == dataset]["cdr3_length"]
                        ax.hist(
                            dataset_lengths,
                            bins=range(int(dataset_lengths.min()), int(dataset_lengths.max()) + 2),
                            alpha=0.6,
                            label=dataset,
                            density=True,
                        )
                    ax.legend(fontsize=8)
                else:
                    # Single distribution
                    ax.hist(
                        chain_data["cdr3_length"],
                        bins=range(int(chain_data["cdr3_length"].min()),
                                 int(chain_data["cdr3_length"].max()) + 2),
                        color="steelblue",
                        alpha=0.7,
                        density=True,
                    )

                ax.set_xlabel("CDR3 Length (nt)", fontsize=11)
                ax.set_ylabel("Density", fontsize=11)
                ax.set_title(f"Chain: {chain}", fontsize=12, fontweight="bold")
                ax.grid(axis="y", alpha=0.3)

            plt.suptitle("CDR3 Length Distribution (Spectratype)",
                        fontsize=14, fontweight="bold", y=1.02)
    else:
        # Single plot without chain separation
        fig, ax = plt.subplots(figsize=figsize)

        if dataset_col and dataset_col in plot_data.columns:
            # Plot by dataset
            for dataset in plot_data[dataset_col].unique():
                dataset_lengths = plot_data[plot_data[dataset_col] == dataset]["cdr3_length"]
                ax.hist(
                    dataset_lengths,
                    bins=range(int(dataset_lengths.min()), int(dataset_lengths.max()) + 2),
                    alpha=0.6,
                    label=dataset,
                    density=True,
                )
            ax.legend(title="Dataset", fontsize=9)
        else:
            # Single distribution
            ax.hist(
                plot_data["cdr3_length"],
                bins=range(int(plot_data["cdr3_length"].min()),
                         int(plot_data["cdr3_length"].max()) + 2),
                color="steelblue",
                alpha=0.7,
                density=True,
            )

        ax.set_xlabel("CDR3 Length (nt)", fontsize=12)
        ax.set_ylabel("Density", fontsize=12)
        ax.set_title("CDR3 Length Distribution (Spectratype)", fontsize=14, fontweight="bold")
        ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()

    if outpath:
        ensure_dir(outpath.parent)
        plt.savefig(outpath, dpi=FIGURE_DPI, bbox_inches="tight")
        logger.info("Saved spectratype plot to %s", outpath)

    plt.close(fig)


def plot_vj_pairing_heatmap(
    receptor_df: pd.DataFrame,
    v_col: str = "v_gene",
    j_col: str = "j_gene",
    top_n: int = 15,
    outpath: Optional[Path] = None,
    figsize: Tuple[int, int] = (10, 8),
) -> None:
    """Generate heatmap of V-J gene pairing frequencies.

    Args:
        receptor_df: DataFrame with receptor data
        v_col: Column name for V gene
        j_col: Column name for J gene
        top_n: Number of top genes to include per segment
        outpath: Output file path for PNG
        figsize: Figure size tuple
    """
    if v_col not in receptor_df.columns or j_col not in receptor_df.columns:
        logger.warning("V or J gene columns not found; skipping V-J pairing plot")
        return

    # Filter valid pairs
    valid_data = receptor_df[[v_col, j_col]].dropna()
    valid_data = valid_data[(valid_data[v_col] != "") & (valid_data[j_col] != "")]

    if valid_data.empty:
        logger.warning("No valid V-J pairs found for heatmap")
        return

    # Get top genes
    top_v = valid_data[v_col].value_counts().head(top_n).index
    top_j = valid_data[j_col].value_counts().head(top_n).index

    # Filter to top genes
    subset = valid_data[valid_data[v_col].isin(top_v) & valid_data[j_col].isin(top_j)]

    # Create contingency table
    pairing_counts = pd.crosstab(subset[v_col], subset[j_col])

    # Sort by frequency
    pairing_counts = pairing_counts.loc[
        pairing_counts.sum(axis=1).sort_values(ascending=False).index,
        pairing_counts.sum(axis=0).sort_values(ascending=False).index,
    ]

    # Create heatmap
    fig, ax = plt.subplots(figsize=figsize)

    sns.heatmap(
        pairing_counts,
        cmap="YlOrRd",
        annot=False,
        fmt="d",
        cbar_kws={"label": "Number of Cells"},
        linewidths=0.5,
        linecolor="gray",
        ax=ax,
    )

    ax.set_xlabel("J Gene", fontsize=12)
    ax.set_ylabel("V Gene", fontsize=12)
    ax.set_title(f"V-J Gene Pairing Frequencies (Top {top_n})", fontsize=14, fontweight="bold")

    # Rotate labels
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=9)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=9)

    plt.tight_layout()

    if outpath:
        ensure_dir(outpath.parent)
        plt.savefig(outpath, dpi=FIGURE_DPI, bbox_inches="tight")
        logger.info("Saved V-J pairing heatmap to %s", outpath)

    plt.close(fig)


def plot_clonal_expansion_summary(
    adata: ad.AnnData,
    clonotype_col: str = "clonotype",
    groupby: Optional[str] = None,
    outpath: Optional[Path] = None,
    figsize: Tuple[int, int] = (10, 6),
) -> None:
    """Generate summary plot of clonal expansion distribution.

    Args:
        adata: AnnData object with clonotype annotations
        clonotype_col: Column name containing clonotype IDs
        groupby: Optional column for grouping (e.g., cancer_type)
        outpath: Output file path for PNG
        figsize: Figure size tuple
    """
    if clonotype_col not in adata.obs.columns:
        logger.warning("Clonotype column not found; skipping expansion summary")
        return

    # Calculate expansion categories
    clonotype_sizes = adata.obs[clonotype_col].map(
        adata.obs[clonotype_col].value_counts()
    ).fillna(1).astype(int)

    expansion_bins = [1, 2, 5, 10, 50, float("inf")]
    expansion_labels = ["1", "2-4", "5-9", "10-49", "≥50"]

    adata.obs["expansion_category"] = pd.cut(
        clonotype_sizes,
        bins=expansion_bins,
        labels=expansion_labels,
        right=False,
    )

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # Plot 1: Overall distribution
    expansion_counts = adata.obs["expansion_category"].value_counts().sort_index()

    axes[0].bar(
        range(len(expansion_counts)),
        expansion_counts.values,
        color="steelblue",
        edgecolor="black",
        alpha=0.7,
    )
    axes[0].set_xticks(range(len(expansion_counts)))
    axes[0].set_xticklabels(expansion_counts.index, fontsize=10)
    axes[0].set_xlabel("Cells per Clonotype", fontsize=11)
    axes[0].set_ylabel("Number of Cells", fontsize=11)
    axes[0].set_title("Overall Clonal Expansion", fontsize=12, fontweight="bold")
    axes[0].grid(axis="y", alpha=0.3)

    # Plot 2: Grouped distribution (if groupby provided)
    if groupby and groupby in adata.obs.columns:
        grouped = adata.obs.groupby([groupby, "expansion_category"]).size().unstack(fill_value=0)
        grouped_norm = grouped.div(grouped.sum(axis=1), axis=0) * 100

        grouped_norm.T.plot(kind="bar", ax=axes[1], colormap="Set3", width=0.8)
        axes[1].set_xlabel("Cells per Clonotype", fontsize=11)
        axes[1].set_ylabel("Percentage of Cells", fontsize=11)
        axes[1].set_title(f"Clonal Expansion by {groupby.replace('_', ' ').title()}",
                         fontsize=12, fontweight="bold")
        axes[1].legend(title=groupby.replace("_", " ").title(),
                      bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=9)
        axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=0)
        axes[1].grid(axis="y", alpha=0.3)
    else:
        axes[1].axis("off")
        axes[1].text(0.5, 0.5, "No grouping variable provided",
                    ha="center", va="center", fontsize=11)

    plt.tight_layout()

    if outpath:
        ensure_dir(outpath.parent)
        plt.savefig(outpath, dpi=FIGURE_DPI, bbox_inches="tight")
        logger.info("Saved clonal expansion summary to %s", outpath)

    plt.close(fig)


def generate_all_tcr_figures(
    adata: ad.AnnData,
    receptor_df: pd.DataFrame,
    diversity_df: Optional[pd.DataFrame],
    figures_dir: Path,
    groupby: str = "cancer_type",
) -> None:
    """Generate complete suite of TCR visualization figures.

    Args:
        adata: AnnData object with clonotype annotations and UMAP
        receptor_df: DataFrame with detailed receptor information
        diversity_df: Optional DataFrame with diversity metrics
        figures_dir: Output directory for figures
        groupby: Grouping variable for stratified analyses
    """
    ensure_dir(figures_dir)
    logger.info("Generating TCR visualization suite in %s", figures_dir)

    # 1. Clonotype frequency barplot
    plot_clonotype_frequency(
        adata,
        clonotype_col="clonotype",
        top_n=20,
        outpath=figures_dir / "clonotype_frequency_top20.png",
        dataset_col="dataset_id" if "dataset_id" in adata.obs.columns else None,
    )

    # 2. Repertoire diversity boxplots
    if diversity_df is not None and not diversity_df.empty:
        plot_repertoire_diversity(
            diversity_df,
            groupby=groupby if groupby in diversity_df.columns else "dataset_id",
            outpath=figures_dir / f"repertoire_diversity_by_{groupby}.png",
        )

    # 3. UMAP colored by clonotype size
    if "X_umap" in adata.obsm:
        plot_umap_clonotype_size(
            adata,
            clonotype_col="clonotype",
            umap_key="X_umap",
            outpath=figures_dir / "umap_clonotype_expansion.png",
        )

    # 4. Spectratype (CDR3 length distribution)
    plot_spectratype(
        receptor_df,
        cdr3_col="cdr3_nt",
        chain_col="chain",
        dataset_col="dataset_id" if "dataset_id" in receptor_df.columns else None,
        outpath=figures_dir / "cdr3_spectratype_by_chain.png",
        separate_chains=True,
    )

    # 5. V-J pairing heatmap
    plot_vj_pairing_heatmap(
        receptor_df,
        v_col="v_gene",
        j_col="j_gene",
        top_n=15,
        outpath=figures_dir / "vj_pairing_heatmap.png",
    )

    # 6. Clonal expansion summary
    plot_clonal_expansion_summary(
        adata,
        clonotype_col="clonotype",
        groupby=groupby if groupby in adata.obs.columns else None,
        outpath=figures_dir / f"clonal_expansion_summary_{groupby}.png",
    )

    logger.info("TCR visualization suite generation complete")
