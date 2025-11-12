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

**Interpretation:**
- **<0.1%**: Normal background (typical for most virome samples)
- **0.1-1%**: Moderate contamination (monitor)
- **>1%**: High contamination (investigate sample)
- **>5%**: Very high contamination (potential sequencing or library prep issue)

#### Visualizations

1. **Bar Plot** (`contamination_bars.png`)
   - Shows contamination levels for each sample
   - Samples sorted by total contamination
   - Helps identify outliers at a glance

2. **Box Plot** (`contamination_boxes.png`)
   - Shows distribution of PhiX and vector contamination across all samples
   - Useful for understanding typical contamination levels

3. **Scatter Plot** (`contamination_scatter.png`)
   - Shows relationship between PhiX and vector contamination
   - Outliers are labeled
   - Can reveal if contamination types co-occur

4. **Heatmap** (`contamination_heatmap.png`)
   - Overview of contamination patterns across all samples
   - Easy to spot samples with elevated contamination

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

### Database

Vector/plasmid database location:
```
workflow/databases/vector_contaminants.fa
```

This should include common contaminants:
- Cloning vectors (pUC, pBR322, etc.)
- Plasmids (pET, pGEX, etc.)
- Synthetic sequences
- Adapter sequences
- UniVec database entries

## Advantages of Flagging vs. Removal

1. **Preserves viral sequences**: Avoids accidentally removing legitimate viral reads that share sequence similarity with contaminants

2. **Better for virome analysis**: VLP-enriched samples should have minimal contamination anyway, so removal is often unnecessary

3. **QC metric**: Contamination levels serve as a quality control indicator:
   - High PhiX: Potential sequencing run issue
   - High vector: Potential library prep contamination
   - Pattern detection: Batch effects, problematic samples

4. **Flexibility**: Researchers can decide post-hoc whether to exclude samples based on contamination levels

5. **Transparency**: All reads retained for downstream analysis, but contamination is documented

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
