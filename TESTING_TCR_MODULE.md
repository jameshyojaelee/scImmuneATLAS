# TCR Module Testing Guide

## Overview

This document describes the comprehensive test suite for the TCR analysis module (`src/atlas/tcr.py`), including fixtures, mocking strategies, and coverage goals.

## Test File Location

```
tests/test_tcr.py          # Main TCR module tests (NEW)
tests/test_tcr_viz.py      # Visualization-specific tests (EXISTING)
tests/test_receptors_io.py # Receptor I/O and utility tests (EXISTING)
```

## Test Coverage Summary

### Covered Components

| Component | Test Functions | Coverage Focus |
|-----------|----------------|----------------|
| **AIRR Schema Conversion** | 5 tests | Schema validation, column mapping, null handling |
| **Scirpy Integration** | 3 tests | Import mocking, merge functionality, error handling |
| **Clonotype Metrics** | 2 tests | Per-dataset metrics, diversity calculations |
| **Repertoire Overlap** | 6 tests | Jaccard index, Morisita-Horn, edge cases |
| **Public Clonotypes** | 4 tests | Detection, thresholds, empty data handling |
| **CDR3 Properties** | 7 tests | Length, charge, hydrophobicity, edge cases |
| **File I/O** | 2 tests | JSON writing, directory creation |
| **Figure Generation** | 3 tests | Matplotlib output, error fallbacks |
| **End-to-End Pipeline** | 2 tests | Full pipeline, disabled configs |
| **Error Handling** | 3 tests | Missing files, empty data, scirpy unavailable |
| **Regression Tests** | 2 tests | Bug prevention, historical issues |

**Total: 39+ test functions**

---

## Fixtures: Synthetic Data Generation

### 1. `mock_receptor_df`

**Purpose**: Generate realistic synthetic TCR receptor dataframe

**Structure**:
```python
@pytest.fixture
def mock_receptor_df():
    """Generate synthetic receptor dataframe with realistic TCR data."""
    np.random.seed(42)
    n_cells = 200

    return pd.DataFrame({
        "cell_id": [f"CELL_{i // 2}" for i in range(n_cells)],
        "dataset_id": np.random.choice(["DATASET_A", "DATASET_B"], n_cells),
        "chain": np.random.choice(["TRA", "TRB"], n_cells),
        "cdr3_nt": [random nucleotide sequences],
        "cdr3_aa": [random amino acid sequences],
        "v_gene": [f"TRBV{i}" for random i],
        "j_gene": [f"TRBJ{i}-{j}" for random i, j],
        "clonotype_id": np.random.choice(clones + singletons),
        # ... other fields
    })
```

**Key Features**:
- 2 chains per cell (TRA/TRB pairs)
- Realistic gene names (TRBV1-30, TRBJ1-2)
- Mix of clonal and singleton cells
- Includes productive/non-productive annotations

---

### 2. `mock_adata_with_tcr`

**Purpose**: Generate AnnData object with scirpy-annotated TCR data

**Structure**:
```python
@pytest.fixture
def mock_adata_with_tcr():
    """Generate synthetic AnnData with TCR annotations from scirpy."""
    adata = ad.AnnData(X=..., obs=pd.DataFrame({
        "cell_type": ["CD4 T", "CD8 T", "NK"],
        "clonotype": [f"clone{i}" for i in range(1, 11)],
        "IR_VJ_1_v_call": [f"TRBV{...}"],
        "IR_VJ_1_j_call": [f"TRBJ{...}"],
        "IR_VJ_1_junction_aa": ["CASS..."],
    }))

    # Add UMAP embedding
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    return adata
```

**Key Features**:
- Includes scirpy-style column names (IR_*)
- UMAP coordinates for visualization tests
- Cell type annotations

---

### 3. `mock_config`

**Purpose**: Generate test configuration dictionary

**Structure**:
```python
@pytest.fixture
def mock_config(tmp_path):
    return {
        "seed": 42,
        "datasets": [
            {
                "id": "DATASET_A",
                "receptor_path": str(tmp_path / "dataset_a_receptors.parquet"),
                "receptor_format": "10x_vdj",
            },
            # ...
        ],
        "tcr": {
            "enabled": True,
            "metrics_dir": str(tmp_path / "metrics"),
            "figures_dir": str(tmp_path / "figures"),
            "min_public_datasets": 2,
        },
        # ...
    }
```

**Key Features**:
- Uses `tmp_path` for isolated file operations
- Includes all required config keys
- Supports per-dataset receptor configs

---

### 4. `mock_datasets_receptor`

**Purpose**: Create dictionary mapping dataset IDs to receptor dataframes

**Structure**:
```python
@pytest.fixture
def mock_datasets_receptor(mock_receptor_df):
    df_a = mock_receptor_df[mock_receptor_df["dataset_id"] == "DATASET_A"]
    df_b = mock_receptor_df[mock_receptor_df["dataset_id"] == "DATASET_B"]

    return {
        "DATASET_A": df_a,
        "DATASET_B": df_b,
    }
```

**Key Features**:
- Splits receptor data by dataset
- Used for multi-dataset overlap analysis

---

## Mocking Scirpy

### Why Mock Scirpy?

Scirpy is an optional dependency that:
1. May not be installed in all test environments
2. Has complex dependencies (e.g., igraph, python-levenshtein)
3. Can be slow for large-scale tests

### Mocking Strategy

#### **Option 1: Skip Tests When Scirpy Unavailable**

```python
@pytest.mark.skipif(
    not hasattr(tcr, "_merge_with_scirpy"),
    reason="Requires scirpy integration"
)
def test_merge_with_scirpy(mock_adata_with_tcr, mock_receptor_df):
    try:
        import scirpy as ir
        scirpy_available = True
    except ImportError:
        scirpy_available = False

    if not scirpy_available:
        pytest.skip("scirpy not available")

    # Test code here...
```

**Pros**:
- Real scirpy behavior when available
- No brittle mocks

**Cons**:
- Tests skipped in minimal environments

---

#### **Option 2: Mock Scirpy Functions Directly**

```python
@patch("src.atlas.tcr._import_scirpy")
def test_tcr_analysis_without_scirpy(mock_import_scirpy, mock_config):
    """Test that TCR analysis fails gracefully when scirpy unavailable."""
    mock_import_scirpy.side_effect = RuntimeError("scirpy is required")

    with pytest.raises(RuntimeError, match="scirpy is required"):
        tcr._import_scirpy()
```

**Pros**:
- Tests run in all environments
- Can test error handling

**Cons**:
- Requires careful mock design

---

#### **Option 3: Mock Scirpy Module at Import**

```python
import sys
from unittest.mock import MagicMock

# Mock entire scirpy module
sys.modules['scirpy'] = MagicMock()

# Now imports of scirpy will use the mock
from src.atlas import tcr
```

**Pros**:
- Prevents import errors
- Full module-level control

**Cons**:
- Can cause confusing behavior if not cleaned up

---

### Recommended Approach

Use **Option 1** (skip tests) for most tests, with **Option 2** (patch specific functions) for error handling tests:

```python
# In conftest.py or test file
try:
    import scirpy as ir
    SCIRPY_AVAILABLE = True
except ImportError:
    SCIRPY_AVAILABLE = False

# Mark tests that require scirpy
pytestmark = pytest.mark.skipif(
    not SCIRPY_AVAILABLE,
    reason="scirpy not installed"
)
```

---

## Test Categories

### 1. AIRR Schema Validation

**Purpose**: Ensure receptor dataframes are correctly converted to AIRR format

**Key Tests**:
- `test_airr_from_dataframe_basic`: Column mapping correctness
- `test_airr_from_dataframe_missing_columns`: Optional field handling
- `test_airr_schema_required_fields`: AIRR compliance
- `test_airr_schema_data_types`: Type checking
- `test_airr_schema_handles_nulls`: Null value handling

**Assertions**:
```python
# Column existence
assert all(col in airr_df.columns for col in required_cols)

# Column mapping
assert airr_df["junction"].equals(mock_receptor_df["cdr3_nt"])
assert airr_df["v_call"].equals(mock_receptor_df["v_gene"])

# Data types
assert airr_df["productive"].dtype == object  # nullable boolean
```

---

### 2. Scirpy Pipeline Functions

**Purpose**: Test integration with scirpy for clonotype definition and metrics

**Key Tests**:
- `test_import_scirpy_missing`: Error when scirpy unavailable
- `test_merge_with_scirpy`: AIRR data merged into AnnData
- `test_compute_dataset_metrics`: Per-dataset metrics computation

**Assertions**:
```python
# Scirpy import
with pytest.raises(RuntimeError, match="scirpy is required"):
    tcr._import_scirpy()

# Metrics structure
assert "n_cells" in metrics
assert "n_clonotypes" in metrics
assert "diversity" in metrics
assert "v_gene_usage" in metrics

# Top clonotypes format
assert isinstance(metrics["top_clonotypes"], list)
assert "clonotype" in metrics["top_clonotypes"][0]
assert "cells" in metrics["top_clonotypes"][0]
```

---

### 3. Repertoire Overlap Analysis

**Purpose**: Test similarity metrics between datasets

**Key Tests**:
- `test_jaccard_index_*`: Jaccard similarity edge cases
- `test_morisita_horn_*`: Morisita-Horn similarity edge cases
- `test_compute_repertoire_overlap`: Pairwise matrix computation
- `test_compute_repertoire_overlap_single_dataset`: Single dataset handling

**Assertions**:
```python
# Jaccard index properties
assert 0.0 <= jaccard <= 1.0
assert jaccard_index(set1, set1) == 1.0  # Self-similarity
assert jaccard_index(disjoint1, disjoint2) == 0.0  # No overlap

# Matrix structure
assert len(overlap_metrics["dataset_ids"]) == n_datasets
assert len(overlap_metrics["jaccard"]) == n_datasets
assert overlap_metrics["jaccard"][i][i] == 1.0  # Diagonal

# Symmetry
assert overlap_metrics["jaccard"][0][1] == overlap_metrics["jaccard"][1][0]
```

---

### 4. Public Clonotype Detection

**Purpose**: Identify clonotypes shared across multiple datasets

**Key Tests**:
- `test_identify_public_clonotypes_basic`: Shared clonotype detection
- `test_identify_public_clonotypes_threshold`: Min dataset threshold
- `test_identify_public_clonotypes_no_sharing`: No shared clonotypes
- `test_identify_public_clonotypes_empty_strings`: Empty ID handling

**Assertions**:
```python
# Structure
assert "n_public_clonotypes" in result
assert "top_public_clonotypes" in result
assert "sharing_distribution" in result

# Public clonotype details
clone_info = next(c for c in result["top_public_clonotypes"] if c["clonotype_id"] == "clone2")
assert clone_info["n_datasets"] == 3
assert set(clone_info["datasets"]) == {"DATASET_A", "DATASET_B", "DATASET_C"}
assert clone_info["total_cells"] == 5

# Threshold behavior
assert result_min2["n_public_clonotypes"] > result_min3["n_public_clonotypes"]
```

---

### 5. CDR3 Property Analysis

**Purpose**: Analyze physicochemical properties of CDR3 sequences

**Key Tests**:
- `test_analyze_cdr3_properties_basic`: Length, charge, hydrophobicity
- `test_analyze_cdr3_properties_varied_lengths`: Length distribution
- `test_analyze_cdr3_properties_charge_calculation`: Charge accuracy
- `test_analyze_cdr3_properties_hydrophobicity`: Hydrophobicity scale
- `test_analyze_cdr3_properties_empty_dataframe`: Empty data handling

**Assertions**:
```python
# Structure
assert "mean_length" in result
assert "median_length" in result
assert "mean_charge" in result
assert "mean_hydrophobicity" in result
assert "length_distribution" in result

# Length calculation
assert result["mean_length"] == 6.0
assert result["median_length"] == 6.0

# Charge calculation (R/K positive, D/E negative)
assert isinstance(result["mean_charge"], float)
assert -1.0 <= result["mean_charge"] <= 1.0  # Normalized

# Empty data fallback
assert result_empty["mean_length"] == 0
assert result_empty["mean_charge"] == 0.0
```

---

### 6. Metrics Output (JSON)

**Purpose**: Test JSON serialization and file writing

**Key Tests**:
- `test_write_json`: Basic JSON writing
- `test_write_json_creates_parent_directories`: Directory creation

**Assertions**:
```python
# File creation
assert output_file.exists()

# Content verification
with open(output_file, "r") as f:
    loaded = json.load(f)
assert loaded == payload
assert loaded["dataset_id"] == "DATASET_A"

# Parent directory creation
nested_path = tmp_path / "nested" / "dir" / "metrics.json"
tcr._write_json(nested_path, payload)
assert nested_path.exists()
```

---

### 7. Figure Generation (Matplotlib)

**Purpose**: Test visualization output using matplotlib testing utilities

**Key Tests**:
- `test_generate_tcr_visualizations`: Basic figure generation
- `test_generate_tcr_visualizations_with_public_clonotypes`: Public clonotype figures
- `test_generate_clonotype_network_fallback`: Error handling

**Assertions**:
```python
# File existence
assert figures_dir.exists()
assert (figures_dir / "repertoire_overlap_jaccard.png").exists()
assert fig_path.stat().st_size > 0  # Non-empty file

# Public clonotype figures (conditional)
if public_clonotypes["n_public_clonotypes"] > 0:
    assert (figures_dir / "public_clonotype_distribution.png").exists()
    assert (figures_dir / "top_public_clonotypes.png").exists()

# Error fallback (should create placeholder)
network_fig = figures_dir / "clonotype_network.png"
if network_fig.exists():
    assert network_fig.stat().st_size > 0
```

**Matplotlib Testing Utilities** (if needed):
```python
import matplotlib.pyplot as plt
from matplotlib.testing.decorators import image_comparison

@image_comparison(baseline_images=['expected_plot'], extensions=['png'])
def test_plot_matches_baseline():
    fig, ax = plt.subplots()
    # Generate plot
    # Pytest will compare against baseline image
```

---

### 8. End-to-End Pipeline Tests

**Purpose**: Test full TCR analysis workflow

**Key Tests**:
- `test_run_tcr_analysis_disabled`: Disabled config handling
- `test_run_tcr_analysis_no_datasets_with_receptors`: No receptor data

**Assertions**:
```python
# Disabled analysis
config["tcr"]["enabled"] = False
result = tcr.run_tcr_analysis(config)
assert result == {}

# No receptor datasets
config["datasets"] = [{"id": "NO_TCR", "path": "..."}]
result = tcr.run_tcr_analysis(config)
assert result == {}
```

---

### 9. Error Handling & Edge Cases

**Purpose**: Test robustness against invalid inputs

**Key Tests**:
- `test_merge_with_scirpy_empty_dataframe`: Empty receptor data
- `test_compute_repertoire_overlap_empty_clonotypes`: Empty clonotype IDs
- `test_load_dataset_artifacts_missing_files`: Missing input files

**Assertions**:
```python
# Empty dataframe handling
adata_merged = tcr._merge_with_scirpy(adata, empty_df)
assert adata_merged.n_obs == adata.n_obs  # Should not crash

# Empty clonotypes
overlap_metrics = _compute_repertoire_overlap(datasets_with_empty_clonotypes)
assert overlap_metrics["jaccard"][0][1] == 0.0

# Missing files
with pytest.raises(FileNotFoundError):
    tcr._load_dataset_artifacts(dataset_cfg, config)
```

---

### 10. Parametrized Tests

**Purpose**: Test functions across multiple parameter values

**Examples**:
```python
@pytest.mark.parametrize("n_datasets", [1, 2, 3, 5])
def test_repertoire_overlap_various_dataset_counts(n_datasets):
    datasets_receptor = {
        f"DATASET_{i}": pd.DataFrame({"clonotype_id": [...]})
        for i in range(n_datasets)
    }
    overlap_metrics = _compute_repertoire_overlap(datasets_receptor)
    assert len(overlap_metrics["dataset_ids"]) == n_datasets


@pytest.mark.parametrize("cdr3_length", [6, 10, 15, 20])
def test_cdr3_properties_various_lengths(cdr3_length):
    receptor_df = pd.DataFrame({"cdr3_aa": ["A" * cdr3_length] * 10})
    result = _analyze_cdr3_properties(receptor_df)
    assert result["mean_length"] == cdr3_length
```

---

## Running Tests

### Basic Test Execution

```bash
# Run all TCR tests
pytest tests/test_tcr.py -v

# Run specific test
pytest tests/test_tcr.py::test_jaccard_index_basic -v

# Run with coverage
pytest tests/test_tcr.py --cov=src.atlas.tcr --cov-report=html
```

### Skip Scirpy-Dependent Tests

```bash
# Skip tests requiring scirpy
pytest tests/test_tcr.py -m "not scirpy_required"
```

### Run Only Fast Tests

```bash
# Skip slow end-to-end tests
pytest tests/test_tcr.py -m "not slow"
```

---

## Coverage Goals

### Target Coverage: 85%+

| Module Component | Target Coverage |
|------------------|----------------|
| AIRR conversion | 95% |
| Scirpy integration | 70% (many branches depend on scirpy availability) |
| Metrics computation | 90% |
| Overlap analysis | 95% |
| Public clonotypes | 95% |
| CDR3 properties | 95% |
| File I/O | 90% |
| Figure generation | 80% (matplotlib outputs hard to test) |
| Error handling | 85% |

### Coverage Report

```bash
pytest tests/test_tcr.py --cov=src.atlas.tcr --cov-report=term-missing

# Generate HTML report
pytest tests/test_tcr.py --cov=src.atlas.tcr --cov-report=html
open htmlcov/index.html
```

---

## Continuous Integration

### GitHub Actions Example

```yaml
name: TCR Module Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: [3.9, 3.10, 3.11]

    steps:
      - uses: actions/checkout@v3

      - name: Set up Python ${{ matrix.python-version }}
        uses: actions/setup-python@v4
        with:
          python-version: ${{ matrix.python-version }}

      - name: Install dependencies
        run: |
          pip install -r requirements.txt
          pip install pytest pytest-cov
          # Optional: Install scirpy
          pip install scirpy || echo "Scirpy not available"

      - name: Run TCR tests
        run: |
          pytest tests/test_tcr.py --cov=src.atlas.tcr --cov-report=xml

      - name: Upload coverage to Codecov
        uses: codecov/codecov-action@v3
        with:
          file: ./coverage.xml
```

---

## Debugging Failed Tests

### Common Issues

#### 1. **Scirpy Import Errors**

**Symptom**: `ModuleNotFoundError: No module named 'scirpy'`

**Solution**:
```bash
pip install scirpy
# Or skip tests
pytest tests/test_tcr.py -k "not scirpy"
```

#### 2. **Random Seed Issues**

**Symptom**: Inconsistent test results

**Solution**:
```python
# Always set seed in fixtures
np.random.seed(42)

# Or use pytest-randomly
pytest tests/test_tcr.py -p no:randomly
```

#### 3. **Matplotlib Backend Issues**

**Symptom**: `RuntimeError: Cannot use Agg backend in non-main thread`

**Solution**:
```python
# In conftest.py
import matplotlib
matplotlib.use('Agg')
```

#### 4. **Temporary Directory Cleanup**

**Symptom**: Disk space issues or stale test files

**Solution**:
```python
# Use tmp_path fixture (auto-cleaned)
def test_example(tmp_path):
    output_dir = tmp_path / "outputs"
    # Will be cleaned up after test
```

---

## Extending the Test Suite

### Adding New Tests

1. **Choose appropriate category** (AIRR validation, metrics, etc.)
2. **Reuse existing fixtures** where possible
3. **Follow assertion patterns** from similar tests
4. **Add parametrize decorators** for robustness
5. **Document expected behavior** in docstring

### Example Template

```python
def test_new_feature(mock_receptor_df, tmp_path):
    """Test description: what does this test verify?

    Expected behavior:
    - Step 1
    - Step 2
    - Edge case handling
    """
    # Arrange
    input_data = mock_receptor_df.copy()

    # Act
    result = tcr.new_function(input_data)

    # Assert
    assert "expected_key" in result
    assert result["expected_value"] > 0
    assert result["list_field"] == expected_list
```

---

## Best Practices

### 1. **Isolation**
- Use `tmp_path` for file I/O
- Don't rely on external files
- Clean up resources (though pytest does this automatically)

### 2. **Determinism**
- Set random seeds (`np.random.seed(42)`)
- Use fixed data for edge cases
- Avoid time-dependent behavior

### 3. **Clarity**
- Descriptive test names (`test_jaccard_index_with_empty_sets`)
- Clear failure messages (`assert result == expected, f"Got {result}, expected {expected}"`)
- Docstrings explaining "why" not "what"

### 4. **Speed**
- Mock expensive operations (scirpy computations)
- Use small datasets (100-200 cells)
- Mark slow tests with `@pytest.mark.slow`

### 5. **Maintainability**
- Reuse fixtures
- Follow repository assertion style
- Keep tests independent (no shared state)

---

## Repository Assertion Style

Based on existing tests in the repository:

```python
# Simple assertions (no message needed)
assert result == expected
assert value > 0
assert "key" in dict

# Existence checks
assert file_path.exists()
assert file_path.stat().st_size > 0

# Type checks
assert isinstance(result, dict)
assert result["field"].dtype == np.float64

# Collection checks
assert len(result) == expected_count
assert set(result) == expected_set
assert result == sorted(result)  # Check ordering

# Regex matches (for errors)
with pytest.raises(ValueError, match="expected pattern"):
    function_that_raises()

# Approximate equality (floats)
assert abs(result - expected) < 1e-6
# Or use pytest.approx
assert result == pytest.approx(expected, rel=1e-5)
```

---

## Summary

- **39+ test functions** covering all major TCR module components
- **Fixtures** generate realistic synthetic TCR data (AIRR-compliant)
- **Scirpy mocking** via skipif decorators and import patches
- **Matplotlib testing** via file existence and size checks
- **Target: 85%+ coverage** with focus on critical paths
- **Repository-style assertions** for consistency

All tests are isolated, deterministic, and follow pytest best practices.
