# TCR Module Tests - Quick Reference

## Running Tests

### Run all TCR tests
```bash
pytest tests/test_tcr.py -v
```

### Run specific test category
```bash
# AIRR schema tests only
pytest tests/test_tcr.py -k "airr" -v

# Overlap analysis tests only
pytest tests/test_tcr.py -k "overlap" -v

# CDR3 property tests only
pytest tests/test_tcr.py -k "cdr3" -v
```

### Run with coverage
```bash
# Terminal report
pytest tests/test_tcr.py --cov=src.atlas.tcr --cov-report=term-missing

# HTML report
pytest tests/test_tcr.py --cov=src.atlas.tcr --cov-report=html
open htmlcov/index.html
```

### Skip scirpy-dependent tests
```bash
pytest tests/test_tcr.py -k "not scirpy"
```

## Test Structure

```
tests/test_tcr.py (NEW - 900+ lines)
├── Fixtures (5 fixtures)
│   ├── mock_receptor_df        # Synthetic TCR receptor data
│   ├── mock_adata_with_tcr     # AnnData with scirpy annotations
│   ├── mock_config              # Test configuration
│   ├── mock_datasets_receptor   # Multi-dataset receptor data
│   └── Uses tmp_path (pytest built-in)
│
├── AIRR Schema Tests (5 tests)
│   ├── test_airr_from_dataframe_basic
│   ├── test_airr_from_dataframe_missing_columns
│   ├── test_airr_schema_required_fields
│   ├── test_airr_schema_data_types
│   └── test_airr_schema_handles_nulls
│
├── Scirpy Integration Tests (3 tests, may skip if unavailable)
│   ├── test_import_scirpy_missing
│   ├── test_merge_with_scirpy
│   └── test_compute_dataset_metrics
│
├── Repertoire Overlap Tests (6 tests)
│   ├── test_jaccard_index_* (4 tests)
│   ├── test_morisita_horn_* (4 tests)
│   ├── test_compute_repertoire_overlap
│   └── test_compute_repertoire_overlap_single_dataset
│
├── Public Clonotype Tests (4 tests)
│   ├── test_identify_public_clonotypes_basic
│   ├── test_identify_public_clonotypes_threshold
│   ├── test_identify_public_clonotypes_no_sharing
│   └── test_identify_public_clonotypes_empty_strings
│
├── CDR3 Property Tests (7 tests)
│   ├── test_analyze_cdr3_properties_basic
│   ├── test_analyze_cdr3_properties_varied_lengths
│   ├── test_analyze_cdr3_properties_charge_calculation
│   ├── test_analyze_cdr3_properties_hydrophobicity
│   ├── test_analyze_cdr3_properties_empty_dataframe
│   └── test_analyze_cdr3_properties_length_distribution
│
├── File I/O Tests (2 tests)
│   ├── test_write_json
│   └── test_write_json_creates_parent_directories
│
├── Figure Generation Tests (3 tests)
│   ├── test_generate_tcr_visualizations
│   ├── test_generate_tcr_visualizations_with_public_clonotypes
│   └── test_generate_clonotype_network_fallback
│
├── End-to-End Tests (2 tests, may skip if scirpy unavailable)
│   ├── test_run_tcr_analysis_disabled
│   └── test_run_tcr_analysis_no_datasets_with_receptors
│
├── Error Handling Tests (3 tests)
│   ├── test_merge_with_scirpy_empty_dataframe
│   ├── test_compute_repertoire_overlap_empty_clonotypes
│   └── test_load_dataset_artifacts_missing_files
│
├── Parametrized Tests (2 tests)
│   ├── test_repertoire_overlap_various_dataset_counts
│   └── test_cdr3_properties_various_lengths
│
└── Regression Tests (2 tests)
    ├── test_clonotype_counts_exclude_empty_strings
    └── test_v_gene_usage_top_n
```

## Expected Output

### Successful run (with scirpy)
```
tests/test_tcr.py::test_airr_from_dataframe_basic PASSED                 [  2%]
tests/test_tcr.py::test_airr_from_dataframe_missing_columns PASSED       [  5%]
tests/test_tcr.py::test_jaccard_index_basic PASSED                       [ 10%]
tests/test_tcr.py::test_morisita_horn_basic PASSED                       [ 15%]
...
==================== 39 passed in 12.45s ====================
```

### Run without scirpy
```
tests/test_tcr.py::test_airr_from_dataframe_basic PASSED                 [  2%]
tests/test_tcr.py::test_merge_with_scirpy SKIPPED (scirpy not available) [  5%]
tests/test_tcr.py::test_compute_dataset_metrics SKIPPED                  [  8%]
...
==================== 34 passed, 5 skipped in 8.32s ====================
```

## Troubleshooting

### Issue: Scirpy tests fail
**Solution**: Install scirpy or skip tests
```bash
pip install scirpy
# OR
pytest tests/test_tcr.py -k "not scirpy"
```

### Issue: Matplotlib backend error
**Solution**: Set backend in conftest.py
```python
import matplotlib
matplotlib.use('Agg')
```

### Issue: Random test failures
**Solution**: Check random seeds are set
```python
np.random.seed(42)  # Should be in fixtures
```

## Adding New Tests

1. Use existing fixtures where possible
2. Follow naming convention: `test_<function>_<scenario>`
3. Add docstring explaining what is tested
4. Use repository assertion style (see TESTING_TCR_MODULE.md)
5. Run with coverage to ensure new code is tested

```python
def test_new_feature(mock_receptor_df, tmp_path):
    """Test new_feature correctly processes receptor data."""
    # Arrange
    input_data = mock_receptor_df.copy()

    # Act
    result = tcr.new_function(input_data)

    # Assert
    assert "key" in result
    assert result["value"] > 0
```

## Coverage Goals

- **Target**: 85%+ overall
- **AIRR conversion**: 95%+
- **Metrics computation**: 90%+
- **Overlap analysis**: 95%+

Run coverage report:
```bash
pytest tests/test_tcr.py --cov=src.atlas.tcr --cov-report=term-missing
```

## Related Test Files

- `tests/test_tcr_viz.py` - Visualization-specific tests (matplotlib output)
- `tests/test_receptors_io.py` - Receptor I/O and utility functions
- `tests/conftest.py` - Shared pytest configuration

## Full Documentation

See `TESTING_TCR_MODULE.md` for comprehensive documentation including:
- Detailed fixture descriptions
- Scirpy mocking strategies
- Assertion patterns
- CI/CD integration examples
