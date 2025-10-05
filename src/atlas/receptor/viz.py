"""Figure generation for immune receptor analytics."""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from .analytics import AnalyticsResult
from .config import global_receptor_config

CLONOTYPE_FREQUENCY_FIG = "clonotype_frequency.png"
UMAP_EXPANSION_FIG = "umap_clonal_expansion.png"
SPECTRATYPE_FIG = "cdr3_spectratype.png"
VJ_USAGE_FIG = "vj_usage_heatmap.png"


sns.set_style("whitegrid")


def _ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def _placeholder(path: Path, message: str) -> None:
    _ensure_parent(path)
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.axis("off")
    ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)
    logging.info("Created placeholder figure at %s", path)


def _plot_clonotype_frequency(result: AnalyticsResult, path: Path) -> None:
    df = result.clonotype_summary
    if df.empty:
        _placeholder(path, "No clonotype assignments available")
        return
    top = df.sort_values("cell_count", ascending=False).head(20)
    top = top.assign(label=top["dataset_id"] + " | " + top["clonotype_id"].astype(str))
    top = top[::-1]  # reverse for horizontal bar order
    fig, ax = plt.subplots(figsize=(8, max(4, len(top) * 0.3)))
    sns.barplot(data=top, x="cell_count", y="label", hue="dataset_id", dodge=False, ax=ax)
    ax.set_xlabel("Cells per clonotype")
    ax.set_ylabel("Dataset | Clonotype")
    ax.legend(title="Dataset", loc="lower right")
    fig.tight_layout()
    _ensure_parent(path)
    fig.savefig(path, dpi=300)
    plt.close(fig)
    logging.info("Saved clonotype frequency plot to %s", path)


def _plot_umap_expansion(config: dict, result: AnalyticsResult, path: Path) -> None:
    atlas_path = Path(global_receptor_config(config).get("atlas_path"))
    if not atlas_path.exists():
        _placeholder(path, "Annotated atlas not found for UMAP overlay")
        return
    adata = ad.read_h5ad(atlas_path)
    if "X_umap" not in adata.obsm:
        _placeholder(path, "UMAP coordinates unavailable")
        return
    embedding = adata.obsm["X_umap"]
    obs = adata.obs
    df = pd.DataFrame(embedding, columns=["UMAP1", "UMAP2"], index=adata.obs_names)
    df["expansion_bin"] = obs.get("clonotype_expansion_bin", pd.Series(index=obs.index)).fillna("unassigned")
    if df.shape[0] > 50000:
        df = df.sample(50000, random_state=0)
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.scatterplot(
        data=df,
        x="UMAP1",
        y="UMAP2",
        hue="expansion_bin",
        s=5,
        linewidth=0,
        palette="viridis",
        ax=ax,
    )
    ax.set_title("Clonal expansion across UMAP")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.legend(title="Expansion", bbox_to_anchor=(1.05, 1), loc="upper left")
    fig.tight_layout()
    _ensure_parent(path)
    fig.savefig(path, dpi=300)
    plt.close(fig)
    logging.info("Saved UMAP expansion overlay to %s", path)


def _plot_spectratype(result: AnalyticsResult, path: Path) -> None:
    df = result.combined
    if df.empty or "cdr3_nt_length" not in df.columns:
        _placeholder(path, "No CDR3 length information available")
        return
    fig, ax = plt.subplots(figsize=(7, 4))
    sns.histplot(data=df, x="cdr3_nt_length", hue="dataset_id", kde=False, element="step", ax=ax, bins=30)
    ax.set_xlabel("CDR3 length (nt)")
    ax.set_ylabel("Cell count")
    ax.set_title("CDR3 spectratype")
    fig.tight_layout()
    _ensure_parent(path)
    fig.savefig(path, dpi=300)
    plt.close(fig)
    logging.info("Saved spectratype plot to %s", path)


def _plot_vj_usage(result: AnalyticsResult, path: Path) -> None:
    df = result.combined
    if df.empty or "v_gene" not in df.columns or "j_gene" not in df.columns:
        _placeholder(path, "No V/J gene usage detected")
        return
    subset = df[(df["v_gene"] != "") & (df["j_gene"] != "")]
    if subset.empty:
        _placeholder(path, "No paired V/J genes available")
        return
    counts = (
        subset.groupby(["v_gene", "j_gene"]).size().reset_index(name="count")
    )
    top_pairs = counts.sort_values("count", ascending=False).head(50)
    pivot = top_pairs.pivot("v_gene", "j_gene", "count").fillna(0)
    fig, ax = plt.subplots(figsize=(max(6, pivot.shape[1] * 0.4), max(5, pivot.shape[0] * 0.3)))
    sns.heatmap(pivot, cmap="mako", ax=ax)
    ax.set_xlabel("J gene")
    ax.set_ylabel("V gene")
    ax.set_title("Top V/J gene usage")
    fig.tight_layout()
    _ensure_parent(path)
    fig.savefig(path, dpi=300)
    plt.close(fig)
    logging.info("Saved V/J usage heatmap to %s", path)


def generate_figures(config: dict, result: AnalyticsResult) -> None:
    figures_dir = result.figures_dir
    _plot_clonotype_frequency(result, figures_dir / CLONOTYPE_FREQUENCY_FIG)
    _plot_umap_expansion(config, result, figures_dir / UMAP_EXPANSION_FIG)
    _plot_spectratype(result, figures_dir / SPECTRATYPE_FIG)
    _plot_vj_usage(result, figures_dir / VJ_USAGE_FIG)
