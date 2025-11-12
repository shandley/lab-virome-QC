# Contamination Flagging System

## Overview

This branch implements a contamination flagging system for the virome QC pipeline. Instead of removing PhiX and vector/plasmid sequences (which could accidentally discard legitimate viral sequences), we now **flag** contamination to generate QC metrics for identifying outlier samples.

## Key Changes

### 1. PhiX Flagging (Replaces Removal)

**Old behavior (`remove_phix` rule):**
- Removed all reads matching PhiX174 reference
- Downstream steps operated on filtered reads
- Risk of removing legitimate viral sequences with similarity to PhiX

**New behavior (`flag_phix` rule):**
- Detects PhiX contamination using BBDuk in stats-only mode
- **Does NOT filter or modify reads**
- Generates statistics file showing % reads matching PhiX
- Located at: `workflow/Snakefile:188-224`

### 2. Vector/Plasmid Flagging (New)

**New rule (`flag_univec`):**
- Checks for common cloning vectors, plasmids, and synthetic sequences
- Uses existing database: `workflow/databases/vector_contaminants.fa`
- Detects contamination without removing reads
- Helps identify laboratory contamination
- Located at: `workflow/Snakefile:226-258`

### 3. Pipeline Flow Updated

**Old flow:**
```
fastp → remove_phix → host_depletion → rrna_removal → viromeqc
```

**New flow:**
```
fastp → [flag_phix + flag_univec] → host_depletion → rrna_removal → viromeqc
         ↓
    contamination_summary + plots
```

- Host depletion now operates directly on `fastp` output
- Contamination flagging runs in parallel (doesn't block main pipeline)
- Added dependency in `host_depletion` to ensure flagging completes first

### 4. New Analysis Scripts

#### `workflow/scripts/aggregate_contamination_stats.py`
- Parses BBDuk stats files for PhiX and vector contamination
- Calculates contamination percentages per sample
- Generates summary table: `results/reports/contamination_summary.tsv`
- Identifies samples with >1% contamination
- Outputs summary statistics (mean, median, max)

#### `workflow/scripts/plot_contamination.py`
- Generates 4 visualization types:
  1. **Bar plot**: Contamination levels per sample (sorted by total)
  2. **Box plot**: Distribution of PhiX and vector contamination
  3. **Scatter plot**: Correlation between PhiX and vector contamination
  4. **Heatmap**: Overview of contamination across all samples
- Highlights outliers with >1% contamination
- Publication-quality figures (300 dpi PNG)

### 5. New Pipeline Outputs

The `rule all` target now includes:
```python
# Contamination flagging summary
f"{OUTDIR}/reports/contamination_summary.tsv",
# Contamination plots
f"{OUTDIR}/reports/contamination_bars.png",
f"{OUTDIR}/reports/contamination_boxes.png",
f"{OUTDIR}/reports/contamination_scatter.png",
f"{OUTDIR}/reports/contamination_heatmap.png"
```

### 6. Read Counting Updated

- Removed `phix_removed` step from `count_reads` rule
- Pipeline now counts: raw → clumpify → fastp → host_depleted → clean
- No intermediate filtering step for PhiX

## Usage

### Running the Pipeline

```bash
# Run the full pipeline (includes contamination flagging)
snakemake --use-conda -j 8

# Run only contamination flagging for existing data
snakemake --use-conda \
    results/reports/contamination_summary.tsv \
    results/reports/contamination_bars.png
```

### Interpreting Results

**Philosophy: Relative Comparison, Not Fixed Thresholds**

The contamination flagging system uses **statistical outlier detection** rather than fixed thresholds. This is important because:
- Normal contamination levels vary between sequencing runs
- PhiX spike-in amounts differ across facilities
- What matters is identifying samples that are *different from the rest of the batch*

**Example:** If 50 samples have 0.5-2% PhiX and one sample has 35% PhiX, that outlier should stand out clearly in the visualizations.

#### Contamination Summary Table

**Location:** `results/reports/contamination_summary.tsv`

**Columns:**
- `sample`: Sample name
- `input_reads`: Total reads analyzed (post-fastp)
- `phix_reads`: Number of reads matching PhiX
- `phix_percent`: % reads matching PhiX
- `vector_reads`: Number of reads matching vectors/plasmids
- `vector_percent`: % reads matching vectors/plasmids
- `total_contamination_percent`: Combined contamination level

**Statistical Outliers (IQR Method):**

The pipeline identifies outliers using the Interquartile Range (IQR) method:
- Calculate Q1 (25th percentile) and Q3 (75th percentile)
- IQR = Q3 - Q1
- Outliers are values > Q3 + 1.5×IQR or < Q1 - 1.5×IQR

This is a standard statistical method that adapts to your data, regardless of whether your run has 0.1% or 5% median contamination.

**Console Output Example:**
```
CONTAMINATION SUMMARY - Identifying Outliers Within This Run
================================================================================
Total samples: 48

PhiX contamination:
  Median:         0.850%
  Mean ± StdDev:  1.203% ± 2.145%
  Range:          0.120% - 12.450%
  25-75th %ile:   0.450% - 1.320%

STATISTICAL OUTLIERS (>1.5 IQR from median):
--------------------------------------------------------------------------------

PhiX outliers (2 samples):
  sample_042                       12.450%  (14.6x median)
  sample_017                        4.230%  (5.0x median)

No vector outliers detected (all samples within expected range)
```

#### Visualizations - Designed for Outlier Detection

All plots use the same statistical outlier detection method and are designed to make deviations from the group visually obvious:

1. **Bar Plot** (`contamination_bars.png`) - **Best for quick visual screening**
   - **Outliers:** Shown in darker colors (orange/magenta vs. yellow/blue)
   - **Outlier labels:** Sample names in red and bold
   - **Reference lines:** Gray dashed line shows median, shaded region shows 25-75th percentile
   - **Use case:** Quick scan to identify which samples deviate from the group

2. **Box Plot** (`contamination_boxes.png`) - **Best for understanding distribution**
   - **Three panels:** PhiX, Vector, and Total contamination
   - **Outliers:** Marked as red diamonds beyond the whiskers
   - **Statistics:** Median and IQR displayed in each panel
   - **Outlier count:** Shown in subplot titles
   - **Use case:** Understand the typical range and identify extreme values

3. **Scatter Plot** (`contamination_scatter.png`) - **Best for correlation analysis**
   - **Outliers:** Shown as red-bordered diamonds with labels
   - **Labels show:** Sample name and fold-change vs. median (e.g., "12.3x")
   - **Reference lines:** Gray dashed lines show median for each axis
   - **Color scale:** Total contamination (darker = higher)
   - **Use case:** See if PhiX and vector contamination correlate, identify samples high in either

4. **Heatmap** (`contamination_heatmap.png`) - **Best for batch overview**
   - **Samples sorted:** Left to right by total contamination (high to low)
   - **Outliers:** Red boxes around outlier columns, red/bold sample names
   - **Values shown:** Actual percentages annotated in each cell
   - **Use case:** Quick overview of the entire batch, spot patterns

## Technical Details

### BBDuk Detection Mode

The flagging rules use BBDuk without output file specification:

```bash
bbduk.sh \
    in1=input_R1.fq in2=input_R2.fq \
    ref=contaminants.fa \
    k=31 hdist=1 \
    stats=output_stats.txt
```

- `k=31`: K-mer size for matching (31bp)
- `hdist=1`: Hamming distance (allows 1 mismatch)
- No `out=` parameter: Reads are not filtered
- `stats=` file contains match counts

### Databases

**Both contamination databases are bundled with the pipeline** - no downloads needed!

**PhiX reference:**
```
workflow/databases/phix174.fasta
```
- PhiX174 complete genome (5,386 bp)
- Standard Illumina sequencing control

**Vector/plasmid database:**
```
workflow/databases/vector_contaminants.fa
```
- Cloning vectors (pUC, pBR322, etc.)
- Plasmids (pET, pGEX, etc.)
- Synthetic sequences
- Adapter sequences
- UniVec database entries

## Advantages of Flagging vs. Removal

1. **Preserves viral sequences**: Avoids accidentally removing legitimate viral reads that share sequence similarity with contaminants

2. **Better for virome analysis**: VLP-enriched samples should have minimal contamination anyway, so removal is often unnecessary

3. **Outlier-based QC**: Statistical outlier detection adapts to your data:
   - No arbitrary fixed thresholds (e.g., ">1% is bad")
   - Identifies samples that deviate from the batch
   - Works regardless of sequencing run contamination baseline
   - Example: In a run where all samples have 5% PhiX, none are flagged as outliers

4. **Interpretable results**: Contamination levels serve as quality indicators:
   - High PhiX outliers: Potential sequencing run issue or incomplete PhiX cleanup
   - High vector outliers: Potential library prep contamination
   - Pattern detection: Batch effects, systematic issues

5. **Flexibility**: Researchers can decide post-hoc whether to exclude samples based on:
   - Statistical outlier status
   - Absolute contamination levels
   - Study-specific tolerance thresholds

6. **Transparency**: All reads retained for downstream analysis, but contamination is documented

## Files Modified

### Core Pipeline
- `workflow/Snakefile`:
  - Replaced `remove_phix` with `flag_phix` (line 188)
  - Added `flag_univec` rule (line 226)
  - Updated `host_depletion` input (line 274)
  - Updated `count_reads` rule (removed phix_removed step)
  - Added `aggregate_contamination_stats` rule (line 470)
  - Added `plot_contamination` rule (line 488)
  - Updated `rule all` targets (line 50)
  - Updated `multiqc` inputs (line 530)

### New Scripts
- `workflow/scripts/aggregate_contamination_stats.py`: Parse and aggregate BBDuk stats
- `workflow/scripts/plot_contamination.py`: Generate contamination visualizations

### Documentation
- `CONTAMINATION_FLAGGING.md`: This file

## Future Enhancements

Possible improvements for future versions:

1. **Thresholds**: Add configurable contamination thresholds to config.yaml
2. **Sample flagging**: Integrate contamination flags into `sample_qc_flags.tsv`
3. **Per-contaminant breakdown**: Show which specific vectors/plasmids matched
4. **MultiQC integration**: Add custom MultiQC module for contamination plots
5. **Trend analysis**: Track contamination across sequencing runs
6. **Additional databases**: Flag adapter dimers, synthetic constructs, etc.

## Questions?

For questions or issues with the contamination flagging system:
1. Check BBDuk stats files: `results/contamination_flagging/*/`
2. Review log files: `results/logs/flag_phix/` and `results/logs/flag_univec/`
3. Verify database exists: `workflow/databases/vector_contaminants.fa`
