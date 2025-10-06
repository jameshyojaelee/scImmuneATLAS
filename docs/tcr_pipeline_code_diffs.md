# TCR Pipeline Integration - Code Diffs

## File 1: Snakefile

### Change 1: Add TCR Rule

**Location**: After `receptor_analytics` rule (line 170)

```diff
if RECEPTOR_ENABLED and RECEPTOR_DATASETS:
    rule receptor_analytics:
        input:
            atlas="processed/integrated_annotated.h5ad",
            qc=[RECEPTOR_QC_MAP[d] for d in RECEPTOR_DATASETS],
            tables=[RECEPTOR_PARQUET_MAP[d] for d in RECEPTOR_DATASETS]
        output:
            summary=str(RECEPTOR_METRICS_DIR / "repertoire_summary.json"),
            expansion=str(RECEPTOR_METRICS_DIR / "clonotype_expansion.tsv"),
            diversity=str(RECEPTOR_METRICS_DIR / "diversity_indices.json"),
            vj=str(RECEPTOR_METRICS_DIR / "vj_usage.tsv"),
            pairing=str(RECEPTOR_METRICS_DIR / "chain_pairing_summary.json"),
            v_usage_json=str(RECEPTOR_METRICS_DIR / "v_usage_summary.json"),
            report_fragment=str(RECEPTOR_METRICS_DIR / "report_section.md"),
            freq_fig=str(RECEPTOR_FIGURES_DIR / "clonotype_frequency.png"),
            umap_fig=str(RECEPTOR_FIGURES_DIR / "umap_clonal_expansion.png"),
            spectratype_fig=str(RECEPTOR_FIGURES_DIR / "cdr3_spectratype.png"),
            vj_fig=str(RECEPTOR_FIGURES_DIR / "vj_usage_heatmap.png"),
        params:
            config_file="config/atlas.yaml"
        shell:
            "python -m src.atlas.receptor --stage analytics --config {params.config_file}"
+
+   rule tcr:
+       input:
+           atlas="processed/integrated_annotated.h5ad",
+           doublet_filtered=expand("data/interim/{dataset}.doublet_filtered.h5ad", dataset=RECEPTOR_DATASETS),
+           tables=[RECEPTOR_PARQUET_MAP[d] for d in RECEPTOR_DATASETS]
+       output:
+           summary=str(RECEPTOR_METRICS_DIR / "tcr_summary.json"),
+           overlap=str(RECEPTOR_METRICS_DIR / "repertoire_overlap.json"),
+           public=str(RECEPTOR_METRICS_DIR / "public_clonotypes.json"),
+           # Main visualization outputs
+           freq_top20=str(RECEPTOR_FIGURES_DIR / "clonotype_frequency_top20.png"),
+           diversity_fig=str(RECEPTOR_FIGURES_DIR / "repertoire_diversity_by_cancer_type.png"),
+           umap_expansion=str(RECEPTOR_FIGURES_DIR / "umap_clonotype_expansion.png"),
+           spectratype_chain=str(RECEPTOR_FIGURES_DIR / "cdr3_spectratype_by_chain.png"),
+           overlap_jaccard=str(RECEPTOR_FIGURES_DIR / "repertoire_overlap_jaccard.png"),
+           vj_heatmap=str(RECEPTOR_FIGURES_DIR / "vj_pairing_heatmap.png"),
+       params:
+           config_file="config/atlas.yaml"
+       shell:
+           "python -m src.atlas.tcr --config {params.config_file}"
```

### Change 2: Update REPORT_INPUTS

**Location**: Before `rule report` (line 235)

```diff
REPORT_INPUTS = [
    "processed/integrated_annotated.h5ad",
    "processed/figures/umap_by_cell_type.png",
]

if RECEPTOR_ENABLED and RECEPTOR_DATASETS:
    REPORT_INPUTS.append(str(RECEPTOR_METRICS_DIR / "repertoire_summary.json"))
    REPORT_INPUTS.append(str(RECEPTOR_METRICS_DIR / "report_section.md"))
+   REPORT_INPUTS.append(str(RECEPTOR_METRICS_DIR / "tcr_summary.json"))
```

## File 2: src/atlas/workflow.py

### Change: Update collect_pipeline_targets()

**Location**: Lines 64-80

```diff
    receptor_cfg = global_receptor_config(config)
    if receptor_cfg.get("enabled"):
        receptor_metrics_dir = Path(
            receptor_cfg.get("qc_metrics_dir", receptor_cfg.get("metrics_dir", "processed/metrics/tcr"))
        )
        receptor_figures_dir = Path(
            receptor_cfg.get("figures_dir", "processed/figures/tcr")
        )
        targets.extend(
            [
                receptor_metrics_dir / "repertoire_summary.json",
                receptor_metrics_dir / "clonotype_expansion.tsv",
                receptor_metrics_dir / "diversity_indices.json",
                receptor_metrics_dir / "vj_usage.tsv",
                receptor_metrics_dir / "chain_pairing_summary.json",
                receptor_metrics_dir / "v_usage_summary.json",
                receptor_metrics_dir / "report_section.md",
                receptor_figures_dir / "clonotype_frequency.png",
                receptor_figures_dir / "umap_clonal_expansion.png",
                receptor_figures_dir / "cdr3_spectratype.png",
                receptor_figures_dir / "vj_usage_heatmap.png",
            ]
        )

-   if receptor_cfg.get("enabled"):
-       tcr_metrics_dir = Path(receptor_cfg.get("qc_metrics_dir", "processed/metrics/tcr"))
-       targets.append(tcr_metrics_dir / "tcr_summary.json")
+   # TCR analysis outputs (scirpy-based)
+   if receptor_cfg.get("enabled"):
+       tcr_metrics_dir = Path(receptor_cfg.get("qc_metrics_dir", "processed/metrics/tcr"))
+       tcr_figures_dir = Path(receptor_cfg.get("figures_dir", "processed/figures/tcr"))
+       targets.extend([
+           # TCR metrics
+           tcr_metrics_dir / "tcr_summary.json",
+           tcr_metrics_dir / "repertoire_overlap.json",
+           tcr_metrics_dir / "public_clonotypes.json",
+           # TCR visualizations
+           tcr_figures_dir / "clonotype_frequency_top20.png",
+           tcr_figures_dir / "repertoire_diversity_by_cancer_type.png",
+           tcr_figures_dir / "umap_clonotype_expansion.png",
+           tcr_figures_dir / "cdr3_spectratype_by_chain.png",
+           tcr_figures_dir / "repertoire_overlap_jaccard.png",
+           tcr_figures_dir / "vj_pairing_heatmap.png",
+       ])
```

## File 3: Makefile

### Change 1: Update .PHONY

**Location**: Line 1

```diff
-.PHONY: setup lint test demo all app env-bootstrap fetch-no-cap run-after-fetch validate-data pipeline receptor
+.PHONY: setup lint test demo all app env-bootstrap fetch-no-cap run-after-fetch validate-data pipeline receptor tcr report
```

### Change 2: Add TCR Target

**Location**: After `receptor` target (lines 31-32)

```diff
receptor:
	scimmuneatlas receptor --config config/atlas.yaml --stage all
+
+tcr:
+	scimmuneatlas tcr --config config/atlas.yaml
```

### Change 3: Update Report Target

**Location**: Line 35

```diff
report:
-	python -m atlas.cli report
+	scimmuneatlas report --config config/atlas.yaml
```

## File 4: src/atlas/cli.py

### No Changes Required

The TCR command was already integrated in previous work:

```python
# Lines 142-155
parser.add_argument("command", choices=[
    "qc",
    "doublets",
    "integrate",
    "annotate",
    "export",
    "viz",
    "benchmark",
    "report",
    "pipeline",
    "validate-data",
    "receptor",
    "tcr",  # ← Already present
])

# Lines 214-215
elif args.command == "tcr":
    tcr.run_tcr_analysis(config)  # ← Already implemented
```

## Summary of Changes

| File | Lines Added | Lines Removed | Net Change |
|------|-------------|---------------|------------|
| Snakefile | 20 | 0 | +20 |
| src/atlas/workflow.py | 17 | 2 | +15 |
| Makefile | 5 | 2 | +3 |
| src/atlas/cli.py | 0 | 0 | 0 |
| **Total** | **42** | **4** | **+38** |

## Verification Commands

### Verify Snakemake Rule Exists

```bash
snakemake --list | grep tcr
# Expected output: tcr
```

### Verify TCR in Pipeline Targets

```python
from src.atlas.workflow import collect_pipeline_targets
from src.atlas.utils import load_config

config = load_config('config/atlas.yaml')
targets = collect_pipeline_targets(config)

tcr_targets = [t for t in targets if 'tcr' in t.lower()]
print(f"Found {len(tcr_targets)} TCR targets:")
for target in tcr_targets:
    print(f"  - {target}")
```

### Verify CLI Command

```bash
scimmuneatlas --help | grep tcr
# Expected output: tcr as one of the commands

scimmuneatlas tcr --help
# Should show usage information
```

### Test Makefile Target

```bash
make -n tcr
# Shows what would be executed (dry-run)

make tcr
# Actually runs TCR analysis
```

## Integration Test

### Full Pipeline with TCR

```bash
# Configure TCR
cat > config/atlas.yaml << 'EOF'
project_name: "Test Atlas"
organism: "human"
seed: 42

tcr:
  enabled: true
  metrics_dir: "processed/metrics/tcr"
  figures_dir: "processed/figures/tcr"

datasets:
  - id: "TEST_01"
    url: "data/raw/test.h5ad"
    platform: "10x"
    cancer_type: "LUAD"
    receptor_path: "data/raw/test_vdj.csv"
    receptor_format: "10x_vdj"

integration:
  method: "scvi"
  batch_key: "dataset_id"
  latent_dim: 30

annotation:
  min_pct: 0.1

outputs:
  metrics_dir: "processed/metrics"
  figures_dir: "processed/figures"
  cellxgene_dir: "processed/cellxgene_release"
EOF

# Run pipeline
snakemake --cores 4

# Verify TCR outputs exist
ls -la processed/metrics/tcr/
ls -la processed/figures/tcr/
```

### Python Fallback Mode

```bash
# Force Python execution
scimmuneatlas pipeline --config config/atlas.yaml --force-python

# Verify same outputs
ls -la processed/metrics/tcr/tcr_summary.json
```

## Rollback Instructions

If you need to revert these changes:

```bash
# Revert Snakefile
git checkout HEAD -- Snakefile

# Revert workflow.py
git checkout HEAD -- src/atlas/workflow.py

# Revert Makefile
git checkout HEAD -- Makefile
```

Or apply these reverse diffs manually.
