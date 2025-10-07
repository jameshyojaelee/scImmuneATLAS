# Streamlit Repertoire Tab - Implementation Guide

## Overview

The Streamlit app has been extended with an optional **Repertoire** tab that provides interactive TCR/BCR analysis when immune receptor data is available.

## Features Implemented

### 1. **Clonotype Distribution Charts**
- **Interactive bar chart** showing top N expanded clonotypes (user-configurable via slider)
- Hover tooltips display clonotype ID and cell count
- Gracefully handles missing data with informative messages

### 2. **Searchable Clonotype Tables**
- **Full clonotype table** with rank, ID, and cell count
- **Real-time search** functionality to filter clonotypes by ID
- **CSV download** button for exporting complete clonotype data
- Displays top 100 by default to avoid performance issues with large datasets

### 3. **V/J Gene Usage Heatmaps**
- **Interactive Plotly heatmap** showing V-J gene pairing frequencies
- User-configurable number of top genes (10-25, default 15)
- Hover tooltips show V gene, J gene, and cell count
- Expandable section with detailed gene usage statistics for top 10 V and J genes

### 4. **Additional Features**
- **Summary metrics**: Total cells with TCR, unique clonotypes, public clonotypes, clonality score
- **CDR3 properties**: Mean/median length, hydrophobicity, length distribution chart
- **Graceful degradation**: Clear messaging when TCR data is unavailable or incomplete

---

## Caching Strategy (`st.cache_data`)

### Why Caching Matters
Loading large single-cell datasets and computing repertoire metrics can be slow. Streamlit's `st.cache_data` decorator ensures expensive operations run only once and are cached for subsequent interactions.

### Functions with Caching

#### 1. **`load_atlas()`** - ✅ Already cached
```python
@st.cache_data
def load_atlas():
    """Load the annotated atlas."""
    adata = ad.read_h5ad(atlas_path)
    return adata
```
- **Purpose**: Avoid re-reading the large H5AD file on every widget interaction
- **Cache key**: Automatically based on file path and modification time

#### 2. **`load_tcr_metrics(config)`** - ✅ New, cached
```python
@st.cache_data
def load_tcr_metrics(config):
    """Load TCR metrics JSON if available."""
    ...
```
- **Purpose**: Load pre-computed TCR summary metrics (clonotype counts, diversity, etc.)
- **Cache invalidation**: Automatically updates when `tcr_summary.json` changes
- **Handles missing files**: Returns `None` gracefully

#### 3. **`load_receptor_data(config)`** - ✅ New, cached
```python
@st.cache_data
def load_receptor_data(config):
    """Load aggregated receptor table (Parquet)."""
    ...
```
- **Purpose**: Load detailed receptor data for V/J heatmaps without re-reading large files
- **Fallback logic**: Searches for aggregated file first, then individual dataset files
- **Handles missing files**: Returns `None` gracefully

### Caching Best Practices Applied

✅ **Immutable return types**: Functions return AnnData, DataFrames, and dicts (Streamlit can hash these)
✅ **No side effects**: Cached functions don't modify global state
✅ **Automatic invalidation**: Cache clears when input files change
✅ **Error handling**: Missing files return `None` instead of crashing, preventing cache pollution

### When NOT to Cache
- UI widget values (sliders, text inputs) - these are already managed by Streamlit session state
- Plotly figure generation - these are fast and depend on user-selected parameters

---

## Handling Missing TCR Data

### Graceful Degradation Strategy

The implementation uses a **three-tier fallback approach**:

#### **Tier 1: Check if TCR data exists in AnnData**
```python
def check_tcr_availability(adata):
    """Check for scirpy-added columns."""
    tcr_columns = ["clonotype", "IR_VJ_1_v_call", "IR_VJ_1_j_call", ...]
    return any(col in adata.obs.columns for col in tcr_columns)
```
- **Outcome**: If no TCR columns found, display warning and hide Repertoire tab

#### **Tier 2: Load pre-computed metrics (optional)**
```python
tcr_metrics = load_tcr_metrics(config)  # Returns None if not found
if tcr_metrics and "global" in tcr_metrics:
    # Display detailed metrics
else:
    st.info("Detailed metrics not available. Showing data from atlas only.")
```
- **Outcome**: Falls back to atlas-only visualization (basic clonotype counts)

#### **Tier 3: Load receptor table for V/J heatmap (optional)**
```python
receptor_df = load_receptor_data(config)  # Returns None if not found
if receptor_df is not None:
    # Show V/J heatmap
else:
    st.info("Receptor data file not found. V/J usage heatmap unavailable.")
```
- **Outcome**: Shows informational message instead of crashing

### User-Facing Messages

**When TCR data is completely missing:**
```
⚠️ TCR/BCR data not available

No immune receptor data was detected in the atlas. This could mean:
- TCR/BCR analysis was not enabled in the pipeline
- The data preprocessing is still in progress
- The integrated atlas was created before receptor analysis

To enable repertoire analysis, ensure TCR/BCR analysis is configured in config/atlas.yaml.
```

**When detailed metrics are missing:**
```
ℹ️ Detailed metrics not available. Showing data from atlas only.
```

**When receptor table is missing:**
```
ℹ️ Receptor data file not found. V/J usage heatmap requires the aggregated receptor table.
```

---

## UI Copywriting Suggestions

### Tab/Section Headers
- **Main header**: "🧬 T-Cell Repertoire Analysis" (clear, concise, emoji for visual interest)
- **Subsections**: Use numbered emojis for visual hierarchy (📊 Repertoire Summary, 🔬 Clonotype Distribution, etc.)

### Metric Cards
```python
st.metric("Total Cells with TCR", "10,234")
st.metric("Unique Clonotypes", "1,523")
st.metric("Public Clonotypes", "42")
st.metric("Clonality", "0.851")  # 1 - (unique_clonotypes / total_cells)
```
- **Rationale**: Large numbers use comma separators for readability
- **Clonality**: Simple metric (0-1 scale) indicating expansion (1 = highly clonal)

### Interactive Widgets
```python
st.slider("Number of top clonotypes to display", 10, 50, 20)
st.text_input("Search clonotypes by ID", "")
st.selectbox("Clonotype column", ["clonotype"])
```
- **Help tooltips** added to selectbox: `help="Select clonotype definition column"`

### Download Buttons
```python
st.download_button(
    label="📥 Download Full Clonotype Table (CSV)",
    data=csv,
    file_name="clonotype_counts.csv",
    mime="text/csv"
)
```
- **Emoji prefix** (📥) makes the button visually distinct
- **Clear file format** in label

### Expandable Sections
```python
with st.expander("📋 Detailed Gene Usage Statistics"):
    # Show V/J gene tables
```
- **Rationale**: Keeps UI clean while providing access to detailed data

### Color Schemes
- **Clonotype bar chart**: Steelblue (professional, consistent with scientific visualizations)
- **V/J heatmap**: YlOrRd (yellow-orange-red, intuitive for frequency data)
- **CDR3 length distribution**: Viridis (perceptually uniform, colorblind-friendly)

---

## Performance Considerations

### Large Dataset Optimization

1. **Display Limits**
   - Show top 100 clonotypes in table by default (full data available via download)
   - Limit V/J heatmap to top 15 genes (user-configurable to 25)

2. **Data Types**
   - Use `pd.read_parquet()` for receptor data (faster than CSV, compressed)
   - Use `st.cache_data` for all file I/O operations

3. **Plot Rendering**
   - Plotly uses WebGL for large scatter plots (better than Matplotlib for interactivity)
   - Heatmaps rendered client-side (no server load after initial data transfer)

### Memory Management
- **AnnData slicing**: Avoid copying entire atlas; use views where possible
- **Garbage collection**: Matplotlib figures explicitly closed with `plt.close(fig)` (not needed here as we use Plotly)

---

## Testing Checklist

### With TCR Data Available
- [ ] Repertoire tab appears in sidebar
- [ ] Summary metrics display correctly
- [ ] Clonotype distribution chart renders
- [ ] Search functionality filters table
- [ ] CSV download works
- [ ] V/J heatmap renders with correct labels
- [ ] CDR3 properties display (if metrics available)

### Without TCR Data
- [ ] Repertoire tab hidden from sidebar
- [ ] Warning message displayed if user navigates directly (shouldn't happen)

### With Partial TCR Data
- [ ] App displays clonotype chart even if metrics JSON missing
- [ ] Graceful fallback message for missing V/J data
- [ ] No error stacktraces visible to user

---

## File Structure

```
app/
└── streamlit_app.py  # Main app with Repertoire tab

processed/
├── integrated_annotated.h5ad           # Main atlas (with TCR columns)
└── aggregated_receptors.parquet        # Optional: Detailed receptor data

metrics/
└── tcr_summary.json                    # Optional: Pre-computed metrics

data/processed/receptors/               # Fallback: Individual receptor files
├── dataset1.parquet
└── dataset2.parquet
```

---

## Future Enhancements

### Potential Additions
1. **Clonotype network visualization** (if scirpy clonotype clustering available)
2. **Public clonotype table** (show clonotypes shared across datasets)
3. **Diversity metric comparison** (Shannon, Simpson, D50 across cancer types)
4. **CDR3 motif search** (find clonotypes matching a user-provided sequence pattern)
5. **Export selected clonotypes** (filter by cancer type, then download)

### API Extensions
If you plan to expose this data via API (e.g., FastAPI):
- `GET /api/clonotypes?top_n=20` - Return top N clonotypes as JSON
- `GET /api/clonotypes/search?query=CASSLQ` - Search by CDR3 sequence
- `GET /api/vj_usage` - Return V/J pairing matrix

---

## Example Usage

### Starting the App
```bash
cd /gpfs/commons/home/jameslee/scImmuneATLAS
streamlit run app/streamlit_app.py
```

### Navigating to Repertoire Tab
1. Load the app (atlas loads automatically)
2. **If TCR data available**: Select "Repertoire" from sidebar dropdown
3. Explore clonotype distribution with interactive slider
4. Search for specific clonotype IDs in the table
5. Download full clonotype data as CSV
6. Examine V/J gene usage heatmap

### Interpreting Metrics
- **Clonality ≈ 0**: Diverse repertoire (many unique clonotypes)
- **Clonality ≈ 1**: Highly clonal (few dominant clonotypes, indicates expansion)
- **Public clonotypes**: Clonotypes shared across multiple samples (potential common antigen response)

---

## Dependencies

All dependencies already present in your environment:
- `streamlit` - Web app framework
- `plotly` - Interactive plots
- `pandas` - Data manipulation
- `anndata` - Single-cell data structure
- `scirpy` - Optional, for TCR analysis (gracefully handles absence)

No additional installations required! 🎉
