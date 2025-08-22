"""Streamlit app for exploring the Single-cell Immune Atlas."""

import logging
from pathlib import Path
from typing import List, Optional

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
import seaborn as sns
import streamlit as st
import yaml

# Configure logging
logging.basicConfig(level=logging.INFO)


@st.cache_data
def load_config():
    """Load configuration file."""
    try:
        with open("config/atlas.yaml", "r") as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        st.error("Config file not found. Please run the pipeline first.")
        st.stop()


@st.cache_data
def load_atlas():
    """Load the annotated atlas."""
    atlas_path = "processed/integrated_annotated.h5ad"
    
    if not Path(atlas_path).exists():
        st.error(f"Atlas file not found at {atlas_path}. Please run the pipeline first.")
        st.stop()
    
    return ad.read_h5ad(atlas_path)


def plot_umap_plotly(adata, color_by, title):
    """Create interactive UMAP plot with Plotly."""
    # Get UMAP coordinates
    umap_coords = adata.obsm["X_umap"]
    
    # Create DataFrame for plotting
    plot_df = pd.DataFrame({
        "UMAP_1": umap_coords[:, 0],
        "UMAP_2": umap_coords[:, 1],
        "Color": adata.obs[color_by].astype(str),
        "Cell_ID": adata.obs_names,
    })
    
    # Add additional metadata for hover
    for col in ["cell_type", "dataset_id", "cancer_type"]:
        if col in adata.obs.columns:
            plot_df[col] = adata.obs[col].astype(str)
    
    # Create plot
    fig = px.scatter(
        plot_df,
        x="UMAP_1",
        y="UMAP_2", 
        color="Color",
        hover_data=["Cell_ID", "cell_type", "dataset_id"],
        title=title,
        width=800,
        height=600
    )
    
    fig.update_traces(marker=dict(size=3, opacity=0.7))
    fig.update_layout(
        xaxis_title="UMAP 1",
        yaxis_title="UMAP 2",
        legend_title=color_by
    )
    
    return fig


def plot_gene_expression(adata, gene, title):
    """Plot gene expression on UMAP."""
    if gene not in adata.var_names:
        st.error(f"Gene {gene} not found in dataset")
        return None
    
    # Get expression values
    gene_idx = adata.var_names.get_loc(gene)
    if hasattr(adata.X, "toarray"):
        expr_values = adata.X[:, gene_idx].toarray().flatten()
    else:
        expr_values = adata.X[:, gene_idx]
    
    # Get UMAP coordinates
    umap_coords = adata.obsm["X_umap"]
    
    # Create DataFrame
    plot_df = pd.DataFrame({
        "UMAP_1": umap_coords[:, 0],
        "UMAP_2": umap_coords[:, 1],
        "Expression": expr_values,
        "Cell_ID": adata.obs_names,
        "Cell_Type": adata.obs["cell_type"].astype(str),
    })
    
    # Create plot
    fig = px.scatter(
        plot_df,
        x="UMAP_1",
        y="UMAP_2",
        color="Expression",
        hover_data=["Cell_ID", "Cell_Type"],
        title=title,
        color_continuous_scale="viridis",
        width=800,
        height=600
    )
    
    fig.update_traces(marker=dict(size=3, opacity=0.7))
    fig.update_layout(
        xaxis_title="UMAP 1",
        yaxis_title="UMAP 2",
        coloraxis_colorbar_title="Expression"
    )
    
    return fig


def plot_proportions(adata, group_by, split_by):
    """Plot cell type proportions."""
    # Create crosstab
    crosstab = pd.crosstab(adata.obs[group_by], adata.obs[split_by])
    
    # Normalize to proportions
    proportions = crosstab.div(crosstab.sum(axis=1), axis=0)
    
    # Create stacked bar plot
    fig = go.Figure()
    
    for cell_type in proportions.columns:
        fig.add_trace(go.Bar(
            name=cell_type,
            x=proportions.index,
            y=proportions[cell_type],
            hovertemplate=f"{cell_type}: %{{y:.1%}}<extra></extra>"
        ))
    
    fig.update_layout(
        barmode="stack",
        title=f"Cell Type Proportions by {group_by}",
        xaxis_title=group_by,
        yaxis_title="Proportion",
        yaxis=dict(tickformat=".0%"),
        width=800,
        height=500
    )
    
    return fig


def show_dataset_summary(adata, config):
    """Show dataset summary statistics."""
    st.subheader("Dataset Summary")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Total Cells", f"{adata.n_obs:,}")
    
    with col2:
        st.metric("Total Genes", f"{adata.n_vars:,}")
    
    with col3:
        st.metric("Cell Types", len(adata.obs["cell_type"].unique()))
    
    # Show dataset composition
    st.subheader("Dataset Composition")
    
    if "dataset_id" in adata.obs.columns:
        dataset_counts = adata.obs["dataset_id"].value_counts()
        st.bar_chart(dataset_counts)
    
    # Show cell type distribution
    st.subheader("Cell Type Distribution")
    cell_type_counts = adata.obs["cell_type"].value_counts()
    st.bar_chart(cell_type_counts)


def main():
    """Main Streamlit app."""
    st.set_page_config(
        page_title="Single-cell Immune Atlas",
        page_icon="🧬",
        layout="wide"
    )
    
    st.title("🧬 Single-cell Immune Atlas Explorer")
    st.markdown("Interactive exploration of tumor-infiltrating lymphocyte scRNA-seq data")
    
    # Load data
    config = load_config()
    adata = load_atlas()
    
    # Sidebar controls
    st.sidebar.header("Visualization Controls")
    
    # Choose visualization type
    viz_type = st.sidebar.selectbox(
        "Visualization Type",
        ["UMAP Overview", "Gene Expression", "Cell Proportions", "Dataset Summary"]
    )
    
    if viz_type == "UMAP Overview":
        st.header("UMAP Overview")
        
        # Choose coloring
        color_options = ["cell_type"]
        if "dataset_id" in adata.obs.columns:
            color_options.append("dataset_id")
        if "cancer_type" in adata.obs.columns:
            color_options.append("cancer_type")
        if "compartment" in adata.obs.columns:
            color_options.append("compartment")
        
        color_by = st.sidebar.selectbox("Color by", color_options)
        
        # Create and show plot
        fig = plot_umap_plotly(adata, color_by, f"UMAP colored by {color_by}")
        st.plotly_chart(fig, use_container_width=True)
        
        # Show color legend info
        st.subheader(f"{color_by} Categories")
        categories = adata.obs[color_by].value_counts()
        st.dataframe(categories.to_frame("Count"))
    
    elif viz_type == "Gene Expression":
        st.header("Gene Expression")
        
        # Gene selection
        gene = st.sidebar.text_input(
            "Gene Symbol", 
            value="CD8A",
            help="Enter a gene symbol (e.g., CD8A, CD4, FOXP3)"
        )
        
        if gene:
            fig = plot_gene_expression(adata, gene, f"{gene} Expression")
            if fig:
                st.plotly_chart(fig, use_container_width=True)
                
                # Show expression statistics
                if gene in adata.var_names:
                    gene_idx = adata.var_names.get_loc(gene)
                    if hasattr(adata.X, "toarray"):
                        expr_values = adata.X[:, gene_idx].toarray().flatten()
                    else:
                        expr_values = adata.X[:, gene_idx]
                    
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Mean Expression", f"{np.mean(expr_values):.3f}")
                    with col2:
                        st.metric("Max Expression", f"{np.max(expr_values):.3f}")
                    with col3:
                        pct_expressing = 100 * np.mean(expr_values > 0)
                        st.metric("% Expressing", f"{pct_expressing:.1f}%")
    
    elif viz_type == "Cell Proportions":
        st.header("Cell Type Proportions")
        
        # Choose grouping variables
        group_options = ["cell_type"]
        if "dataset_id" in adata.obs.columns:
            group_options.append("dataset_id")
        if "cancer_type" in adata.obs.columns:
            group_options.append("cancer_type")
        
        if len(group_options) >= 2:
            group_by = st.sidebar.selectbox("Group by", group_options[1:])
            split_by = st.sidebar.selectbox("Split by", ["cell_type"])
            
            fig = plot_proportions(adata, group_by, split_by)
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("Need at least 2 grouping variables for proportion analysis")
    
    elif viz_type == "Dataset Summary":
        show_dataset_summary(adata, config)
    
    # Footer
    st.sidebar.markdown("---")
    st.sidebar.markdown("**Single-cell Immune Atlas**")
    st.sidebar.markdown(f"Version: {config.get('version', '0.1.0')}")


if __name__ == "__main__":
    main()
