"""Streamlit app for exploring the Single-cell Immune Atlas."""

# Ensure we do not import user site-packages (avoids NumPy/Matplotlib mismatches)
import os
import sys
try:
    import site  # type: ignore
    user_site = site.getusersitepackages() if hasattr(site, "getusersitepackages") else None
    if user_site and user_site in sys.path:
        sys.path.remove(user_site)
except Exception:
    pass
os.environ.setdefault("PYTHONNOUSERSITE", "1")

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

try:
    import scirpy as ir
    SCIRPY_AVAILABLE = True
except ImportError:
    SCIRPY_AVAILABLE = False

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

    adata = ad.read_h5ad(atlas_path)
    # Convenience: create a symbol column for gene lookup if possible
    if "feature_name" in adata.var.columns and "symbol" not in adata.var.columns:
        adata.var["symbol"] = adata.var["feature_name"].astype(str)
    return adata


@st.cache_data
def load_tcr_metrics(config):
    """Load TCR metrics if available."""
    metrics_dir = Path(config.get("outputs", {}).get("metrics_dir", "metrics"))
    summary_file = metrics_dir / "tcr_summary.json"

    if summary_file.exists():
        import json
        with open(summary_file, "r") as f:
            return json.load(f)
    return None


@st.cache_data
def load_receptor_data(config):
    """Load aggregated receptor table if available."""
    processed_dir = Path(config.get("outputs", {}).get("processed_dir", "processed"))
    receptor_file = processed_dir / "aggregated_receptors.parquet"

    if receptor_file.exists():
        return pd.read_parquet(receptor_file)

    # Fallback: try to find individual dataset receptor files
    receptor_dir = Path("data/processed/receptors")
    if receptor_dir.exists():
        parquet_files = list(receptor_dir.glob("*.parquet"))
        if parquet_files:
            dfs = []
            for f in parquet_files:
                df = pd.read_parquet(f)
                df["dataset_id"] = f.stem
                dfs.append(df)
            return pd.concat(dfs, ignore_index=True)

    return None


def check_tcr_availability(adata):
    """Check if TCR/BCR data is available in the AnnData object."""
    # Check for scirpy-added columns
    tcr_columns = ["clonotype", "IR_VJ_1_v_call", "IR_VJ_1_j_call", "IR_VJ_1_junction_aa"]
    has_tcr = any(col in adata.obs.columns for col in tcr_columns)
    return has_tcr


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


def plot_clonotype_distribution(adata, clonotype_col="clonotype", top_n=20):
    """Plot top clonotype distribution with Plotly."""
    if clonotype_col not in adata.obs.columns:
        st.warning(f"Clonotype column '{clonotype_col}' not found in data")
        return None

    # Get clonotype counts
    clonotype_counts = adata.obs[clonotype_col].value_counts().head(top_n)

    if clonotype_counts.empty:
        st.info("No clonotypes detected")
        return None

    # Create bar chart
    fig = go.Figure(data=[
        go.Bar(
            x=list(range(1, len(clonotype_counts) + 1)),
            y=clonotype_counts.values,
            text=clonotype_counts.values,
            textposition="auto",
            hovertemplate="<b>Rank %{x}</b><br>Clonotype: %{customdata}<br>Cells: %{y}<extra></extra>",
            customdata=clonotype_counts.index,
            marker=dict(color="steelblue")
        )
    ])

    fig.update_layout(
        title=f"Top {top_n} Expanded Clonotypes",
        xaxis_title="Clonotype Rank",
        yaxis_title="Number of Cells",
        height=500,
        showlegend=False
    )

    return fig


def plot_vj_usage_heatmap(receptor_df, v_col="v_gene", j_col="j_gene", top_n=15):
    """Create V/J gene usage heatmap with Plotly."""
    if receptor_df is None or receptor_df.empty:
        st.warning("No receptor data available for V/J usage analysis")
        return None

    if v_col not in receptor_df.columns or j_col not in receptor_df.columns:
        st.warning(f"V or J gene columns not found in receptor data")
        return None

    # Filter valid pairs
    valid_data = receptor_df[[v_col, j_col]].dropna()
    valid_data = valid_data[(valid_data[v_col] != "") & (valid_data[j_col] != "")]

    if valid_data.empty:
        st.info("No valid V-J pairs found")
        return None

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
    fig = go.Figure(data=go.Heatmap(
        z=pairing_counts.values,
        x=pairing_counts.columns,
        y=pairing_counts.index,
        colorscale="YlOrRd",
        hovertemplate="V: %{y}<br>J: %{x}<br>Cells: %{z}<extra></extra>",
        colorbar=dict(title="Cells")
    ))

    fig.update_layout(
        title=f"V-J Gene Pairing Frequencies (Top {top_n})",
        xaxis_title="J Gene",
        yaxis_title="V Gene",
        height=600,
        xaxis=dict(tickangle=45),
    )

    return fig


def show_repertoire_tab(adata, config):
    """Display the Repertoire analysis tab."""
    st.header("🧬 T-Cell Repertoire Analysis")

    # Check if TCR data is available
    has_tcr = check_tcr_availability(adata)

    if not has_tcr:
        st.warning("""
        ⚠️ **TCR/BCR data not available**

        No immune receptor data was detected in the atlas. This could mean:
        - TCR/BCR analysis was not enabled in the pipeline
        - The data preprocessing is still in progress
        - The integrated atlas was created before receptor analysis

        To enable repertoire analysis, ensure TCR/BCR analysis is configured in `config/atlas.yaml`.
        """)
        return

    # Load additional TCR metrics
    tcr_metrics = load_tcr_metrics(config)
    receptor_df = load_receptor_data(config)

    # Summary metrics
    st.subheader("📊 Repertoire Summary")

    if tcr_metrics and "global" in tcr_metrics:
        global_metrics = tcr_metrics["global"]
        col1, col2, col3, col4 = st.columns(4)

        with col1:
            st.metric("Total Cells with TCR", f"{global_metrics.get('total_cells', 0):,}")

        with col2:
            st.metric("Unique Clonotypes", f"{global_metrics.get('total_clonotypes', 0):,}")

        with col3:
            if "n_public_clonotypes" in global_metrics:
                st.metric("Public Clonotypes", f"{global_metrics['n_public_clonotypes']:,}")

        with col4:
            # Calculate clonality (simple metric)
            if global_metrics.get('total_clonotypes', 0) > 0:
                clonality = 1 - (global_metrics['total_clonotypes'] / global_metrics.get('total_cells', 1))
                st.metric("Clonality", f"{clonality:.3f}")
    else:
        st.info("Detailed metrics not available. Showing data from atlas only.")

    # Clonotype distribution
    st.subheader("🔬 Clonotype Distribution")

    clonotype_col = st.selectbox(
        "Clonotype column",
        ["clonotype"] if "clonotype" in adata.obs.columns else [],
        help="Select clonotype definition column"
    )

    if clonotype_col:
        top_n = st.slider("Number of top clonotypes to display", 10, 50, 20)
        fig = plot_clonotype_distribution(adata, clonotype_col, top_n)
        if fig:
            st.plotly_chart(fig, use_container_width=True)

        # Searchable clonotype table
        st.subheader("🔍 Top Clonotypes (Searchable)")

        clonotype_counts = adata.obs[clonotype_col].value_counts().reset_index()
        clonotype_counts.columns = ["Clonotype ID", "Cell Count"]
        clonotype_counts["Rank"] = range(1, len(clonotype_counts) + 1)
        clonotype_counts = clonotype_counts[["Rank", "Clonotype ID", "Cell Count"]]

        # Add search functionality
        search_term = st.text_input("Search clonotypes by ID", "")
        if search_term:
            filtered_table = clonotype_counts[
                clonotype_counts["Clonotype ID"].astype(str).str.contains(search_term, case=False, na=False)
            ]
            st.dataframe(filtered_table, use_container_width=True)
        else:
            st.dataframe(clonotype_counts.head(100), use_container_width=True)

        # Download button
        csv = clonotype_counts.to_csv(index=False)
        st.download_button(
            label="📥 Download Full Clonotype Table (CSV)",
            data=csv,
            file_name="clonotype_counts.csv",
            mime="text/csv"
        )

    # V/J usage heatmap
    st.subheader("📈 V-J Gene Usage Heatmap")

    if receptor_df is not None:
        top_n_vj = st.slider("Number of top V/J genes to display", 10, 25, 15)
        fig_vj = plot_vj_usage_heatmap(receptor_df, top_n=top_n_vj)
        if fig_vj:
            st.plotly_chart(fig_vj, use_container_width=True)

        # Additional gene usage stats
        with st.expander("📋 Detailed Gene Usage Statistics"):
            col1, col2 = st.columns(2)

            with col1:
                st.markdown("**Top 10 V Genes**")
                v_usage = receptor_df["v_gene"].value_counts().head(10)
                st.dataframe(v_usage.to_frame("Count"))

            with col2:
                st.markdown("**Top 10 J Genes**")
                j_usage = receptor_df["j_gene"].value_counts().head(10)
                st.dataframe(j_usage.to_frame("Count"))
    else:
        st.info("Receptor data file not found. V/J usage heatmap requires the aggregated receptor table.")

    # CDR3 properties (if available in metrics)
    if tcr_metrics and "global" in tcr_metrics and "cdr3_properties" in tcr_metrics["global"]:
        st.subheader("🧪 CDR3 Sequence Properties")

        cdr3_props = tcr_metrics["global"]["cdr3_properties"]

        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Mean CDR3 Length", f"{cdr3_props.get('mean_length', 0):.1f} AA")
        with col2:
            st.metric("Median CDR3 Length", f"{cdr3_props.get('median_length', 0):.1f} AA")
        with col3:
            st.metric("Mean Hydrophobicity", f"{cdr3_props.get('mean_hydrophobicity', 0):.3f}")

        # Length distribution
        if "length_distribution" in cdr3_props:
            length_dist = pd.DataFrame(
                list(cdr3_props["length_distribution"].items()),
                columns=["CDR3 Length (AA)", "Count"]
            ).sort_values("CDR3 Length (AA)")

            fig_length = px.bar(
                length_dist,
                x="CDR3 Length (AA)",
                y="Count",
                title="CDR3 Length Distribution",
                labels={"Count": "Number of Sequences"}
            )
            st.plotly_chart(fig_length, use_container_width=True)


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
    viz_options = ["UMAP Overview", "Gene Expression", "Cell Proportions", "Dataset Summary"]

    # Add Repertoire tab if TCR data is available
    if check_tcr_availability(adata):
        viz_options.append("Repertoire")

    viz_type = st.sidebar.selectbox(
        "Visualization Type",
        viz_options
    )
    
    if viz_type == "UMAP Overview":
        st.header("UMAP Overview")
        
        # Choose coloring
        color_options = ["cell_type"]
        if "cell_type_pred" in adata.obs.columns:
            color_options.append("cell_type_pred")
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
        
        # Gene selection with symbol fallback
        gene = st.sidebar.text_input(
            "Gene Symbol or ID", 
            value="CD8A",
            help="Enter a gene symbol (e.g., CD8A) or an ID present in var_names"
        )
        
        def resolve_gene_key(g: str) -> Optional[str]:
            if g in adata.var_names:
                return g
            if "symbol" in adata.var.columns:
                matches = adata.var.index[adata.var["symbol"].astype(str) == g]
                if len(matches) > 0:
                    return matches[0]
            return None
        
        if gene:
            key = resolve_gene_key(gene)
            if key is None:
                st.error(f"Gene {gene} not found in var_names or symbol column")
            else:
                fig = plot_gene_expression(adata, key, f"{gene} Expression")
                if fig:
                    st.plotly_chart(fig, use_container_width=True)
                
                # Show expression statistics
                if key in adata.var_names:
                    gene_idx = adata.var_names.get_loc(key)
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

    elif viz_type == "Repertoire":
        show_repertoire_tab(adata, config)

    # Footer
    st.sidebar.markdown("---")
    st.sidebar.markdown("**Single-cell Immune Atlas**")
    st.sidebar.markdown(f"Version: {config.get('version', '0.1.0')}")


if __name__ == "__main__":
    main()
