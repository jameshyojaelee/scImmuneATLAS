# TCR/BCR Report Integration Summary

## Overview

Modified `src/atlas/utils.py::generate_report()` to automatically insert a comprehensive "TCR/BCR Repertoire Analysis" section with metrics tables and figure thumbnails. The implementation degrades gracefully when TCR is disabled or data is unavailable.

## Changes Made

### 1. New Function: `_generate_tcr_section()`

**Location**: `src/atlas/utils.py` (lines 252-470)

**Purpose**: Generate TCR/BCR section for the markdown report

**Key Features**:
- Checks if TCR is enabled in config
- Loads metrics from `processed/metrics/tcr/tcr_summary.json`
- Creates formatted markdown tables with repertoire statistics
- Embeds figure thumbnails with proper relative paths
- Returns empty list if disabled/unavailable (graceful degradation)

**Function Signature**:
```python
def _generate_tcr_section(
    config: Dict,
    metrics_dir: Path,
    figures_dir: Path,
    embed_figures: bool
) -> List[str]:
    """Generate TCR/BCR repertoire section for report.

    Returns a list of markdown lines. Returns empty list if
    TCR is disabled or no data available.
    """
```

### 2. Modified Function: `generate_report()`

**Location**: `src/atlas/utils.py` (lines 588-591)

**Change**: Replaced old fragment-based approach with new comprehensive section

**Before**:
```python
receptor_cfg = config.get("tcr") or config.get("receptor", {}) or {}
fragment_dir = Path(
    receptor_cfg.get(
        "qc_metrics_dir",
        receptor_cfg.get("metrics_dir", "processed/metrics/tcr"),
    )
)
fragment_path = fragment_dir / "report_section.md"
if fragment_path.exists():
    fragment_text = fragment_path.read_text().strip()
    if fragment_text:
        lines.append(fragment_text)
        lines.append("")
```

**After**:
```python
# Generate TCR/BCR section (replaces old fragment-based approach)
tcr_section = _generate_tcr_section(config, metrics_dir, figures_dir, embed_figures)
if tcr_section:
    lines.extend(tcr_section)
```

## Report Section Structure

The generated TCR section includes:

### 1. Global Repertoire Statistics (Table)
- Total cells with receptor data
- Unique clonotypes
- Public clonotypes count
- Mean repertoire overlap (Jaccard index)
- Mean repertoire similarity (Morisita-Horn index)
- CDR3 properties (length, charge, hydrophobicity)

### 2. Per-Dataset Repertoire Metrics (Table)
Columns:
- Dataset ID
- Number of cells
- Number of clonotypes
- Shannon entropy
- D50 (diversity index)
- Top clonotype size

### 3. V/J Gene Usage (Two Tables)
- Top 10 V genes with cell counts
- Top 10 J genes with cell counts

### 4. Public Clonotypes (Table)
- Summary count
- Top 5 public clonotypes with:
  - Clonotype ID
  - Datasets where found
  - Total cell count

### 5. Repertoire Visualizations (Embedded Images)
If `embed_figures: true` in config:
- Top 20 expanded clonotypes
- Repertoire diversity by cancer type
- UMAP colored by clonotype size
- CDR3 length distribution (spectratype)
- Repertoire overlap (Jaccard index)
- V-J gene pairing frequencies

### 6. Additional Metrics (File References)
- Links to JSON metric files
- Per-dataset metric file count

## Data Sources

The function reads from these files (all paths configurable):

```
processed/metrics/tcr/
├── tcr_summary.json           # Main metrics file (required)
├── repertoire_overlap.json    # Overlap matrices (optional)
├── public_clonotypes.json     # Public clonotype data (optional)
├── DATASET_A_tcr_metrics.json # Per-dataset metrics (optional)
├── DATASET_B_tcr_metrics.json
└── ...

processed/figures/tcr/
├── clonotype_frequency_top20.png
├── repertoire_diversity_by_cancer_type.png
├── umap_clonotype_expansion.png
├── cdr3_spectratype_by_chain.png
├── repertoire_overlap_jaccard.png
└── vj_pairing_heatmap.png
```

## Graceful Degradation

### Scenario 1: TCR Disabled
```yaml
# config/atlas.yaml
tcr:
  enabled: false
```
**Result**: TCR section completely omitted from report.

### Scenario 2: TCR Enabled, No Data
```python
# tcr_summary.json doesn't exist
```
**Result**: TCR section completely omitted (analysis hasn't run yet).

### Scenario 3: Partial Metrics
```python
# Some metrics missing from tcr_summary.json
{
  "global": {
    "total_cells": 1000,
    # "total_clonotypes" missing
  }
}
```
**Result**: Table shows available metrics, omits missing ones.

### Scenario 4: No Figures
```python
# figures/tcr/ directory empty
```
**Result**:
```markdown
### Repertoire Visualizations

_TCR/BCR figures not yet generated._
```

### Scenario 5: JSON Error
```python
# Corrupted JSON file
```
**Result**:
```markdown
## TCR/BCR Repertoire Analysis

_TCR/BCR analysis metrics could not be loaded._
_Error: JSONDecodeError: Expecting value: line 1 column 1 (char 0)_
```

## Configuration Example

```yaml
# config/atlas.yaml
tcr:
  enabled: true
  metrics_dir: "processed/metrics/tcr"
  figures_dir: "processed/figures/tcr"

report:
  embed_figures: true  # Set to false to omit figures
```

## Example Usage

### Running Report Generation

```python
from atlas.utils import generate_report

# After running full pipeline including TCR analysis
generate_report("config/atlas.yaml")

# Output: processed/report.md
```

### Command Line

```bash
# Generate report after pipeline completion
python -m atlas.utils --config config/atlas.yaml

# Or via the workflow
atlas-cli pipeline --config config/atlas.yaml
```

## Code Snippets

### Loading TCR Summary

```python
tcr_summary_path = tcr_metrics_dir / "tcr_summary.json"
if tcr_summary_path.exists():
    tcr_summary = json.loads(tcr_summary_path.read_text())
    global_metrics = tcr_summary.get("global", {})
```

### Creating Markdown Tables

```python
lines.append("| Metric | Value |")
lines.append("|---|---|")
lines.append(f"| Total cells | {global_metrics['total_cells']:,} |")
```

### Embedding Figures

```python
for fig_name, title in tcr_figure_names:
    fig_path = tcr_figures_dir / fig_name
    if fig_path.exists():
        rel_path = os.path.relpath(fig_path, report_dir)
        lines.append(f"#### {title}")
        lines.append(f"![{title}]({rel_path})")
```

### Handling Missing Data

```python
shannon = diversity.get("shannon_entropy", "N/A")
shannon_str = f"{shannon:.2f}" if isinstance(shannon, (int, float)) else shannon
```

## Testing

### Manual Test

```python
# Create mock metrics
import json
from pathlib import Path

metrics_dir = Path("processed/metrics/tcr")
metrics_dir.mkdir(parents=True, exist_ok=True)

mock_summary = {
    "global": {
        "total_cells": 1000,
        "total_clonotypes": 500,
    },
    "datasets": {
        "TEST_DS": {
            "n_cells": 1000,
            "n_clonotypes": 500,
            "diversity": {"shannon_entropy": 5.5},
        }
    }
}

(metrics_dir / "tcr_summary.json").write_text(json.dumps(mock_summary))

# Generate report
from atlas.utils import generate_report
generate_report("config/atlas.yaml")

# Check output
print(Path("processed/report.md").read_text())
```

## Benefits

1. **Automated**: No manual report editing required
2. **Comprehensive**: Covers all TCR metrics in one section
3. **Visual**: Embeds publication-quality figures
4. **Robust**: Handles missing data gracefully
5. **Configurable**: Respects config settings
6. **Maintainable**: Centralized in one function
7. **Consistent**: Follows existing report formatting

## File Statistics

- Modified: `src/atlas/utils.py` (+220 lines)
- Created: `docs/example_tcr_report_section.md` (example output)
- Created: `docs/tcr_report_integration_summary.md` (this file)
