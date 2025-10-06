# Example TCR/BCR Report Section

This is an example of what the automatically generated TCR/BCR section looks like in the atlas report (`processed/report.md`).

---

## TCR/BCR Repertoire Analysis

### Global Repertoire Statistics

| Metric | Value |
|---|---|
| Total cells with receptor data | 12,458 |
| Unique clonotypes | 8,234 |
| Public clonotypes (shared across datasets) | 156 |
| Mean repertoire overlap (Jaccard) | 0.127 |
| Mean repertoire similarity (Morisita-Horn) | 0.342 |
| Mean CDR3 length (AA) | 14.2 |
| Mean CDR3 charge | 0.35 |
| Mean CDR3 hydrophobicity | -0.82 |

### Per-Dataset Repertoire Metrics

| Dataset | Cells | Clonotypes | Shannon Entropy | D50 | Top Clonotype Size |
|---|---|---|---|---|---|
| BRCA_001 | 3,421 | 2,156 | 6.84 | 89 | 245 |
| COAD_002 | 2,987 | 2,034 | 7.12 | 112 | 189 |
| LUAD_003 | 4,123 | 3,012 | 7.45 | 134 | 167 |
| STAD_004 | 1,927 | 1,032 | 5.98 | 67 | 312 |

### V/J Gene Usage

Top 10 most frequent V genes:

| V Gene | Cells |
|---|---|
| TRBV7-2 | 1,234 |
| TRBV20-1 | 987 |
| TRBV5-1 | 856 |
| TRBV12-3 | 743 |
| TRBV19 | 698 |
| TRBV6-5 | 654 |
| TRBV28 | 612 |
| TRBV9 | 589 |
| TRBV7-9 | 534 |
| TRBV29-1 | 487 |

Top 10 most frequent J genes:

| J Gene | Cells |
|---|---|
| TRBJ2-7 | 2,456 |
| TRBJ2-1 | 2,123 |
| TRBJ1-1 | 1,987 |
| TRBJ2-3 | 1,654 |
| TRBJ1-2 | 1,432 |
| TRBJ2-5 | 1,289 |
| TRBJ1-4 | 987 |
| TRBJ2-2 | 876 |
| TRBJ1-6 | 765 |
| TRBJ2-6 | 543 |

### Public Clonotypes

Identified **156** clonotypes shared across multiple datasets.

Top 5 most prevalent public clonotypes:

| Clonotype ID | Datasets | Total Cells |
|---|---|
| `clonotype_00342` | BRCA_001, COAD_002, LUAD_003, STAD_004 | 456 |
| `clonotype_01289` | BRCA_001, LUAD_003, STAD_004 | 287 |
| `clonotype_00567` | COAD_002, LUAD_003, STAD_004 | 234 |
| `clonotype_02134` | BRCA_001, COAD_002, LUAD_003 | 189 |
| `clonotype_00891` | BRCA_001, STAD_004 | 156 |

### Repertoire Visualizations

#### Top 20 Expanded Clonotypes
![Top 20 Expanded Clonotypes](figures/tcr/clonotype_frequency_top20.png)

#### Repertoire Diversity by Cancer Type
![Repertoire Diversity by Cancer Type](figures/tcr/repertoire_diversity_by_cancer_type.png)

#### UMAP Colored by Clonotype Size
![UMAP Colored by Clonotype Size](figures/tcr/umap_clonotype_expansion.png)

#### CDR3 Length Distribution (Spectratype)
![CDR3 Length Distribution (Spectratype)](figures/tcr/cdr3_spectratype_by_chain.png)

#### Repertoire Overlap (Jaccard Index)
![Repertoire Overlap (Jaccard Index)](figures/tcr/repertoire_overlap_jaccard.png)

#### V-J Gene Pairing Frequencies
![V-J Gene Pairing Frequencies](figures/tcr/vj_pairing_heatmap.png)

### Additional TCR/BCR Metrics

- **Summary metrics**: `processed/metrics/tcr/tcr_summary.json`
- **Repertoire overlap**: `processed/metrics/tcr/repertoire_overlap.json`
- **Public clonotypes**: `processed/metrics/tcr/public_clonotypes.json`
- **Per-dataset metrics**: 4 files in `processed/metrics/tcr/`

---

## Notes on Graceful Degradation

The report generation handles missing data gracefully:

1. **TCR disabled**: If `tcr.enabled: false` in config, the entire section is omitted.

2. **No TCR data**: If TCR is enabled but `tcr_summary.json` doesn't exist, the section is omitted.

3. **Partial data**: If some metrics are missing:
   ```markdown
   | Shannon Entropy | N/A |
   ```

4. **No figures**: If figures haven't been generated:
   ```markdown
   _TCR/BCR figures not yet generated._
   ```

5. **JSON parsing errors**: If metrics files are corrupted:
   ```markdown
   _TCR/BCR analysis metrics could not be loaded._
   _Error: JSONDecodeError: Expecting value: line 1 column 1 (char 0)_
   ```
