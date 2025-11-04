# Data Directory

This directory is for storing your raw sequencing data and intermediate files during pipeline execution.

## Directory Structure

```
data/
├── raw/              # Place your raw FASTQ files here
│   ├── sample1_R1.fastq.gz
│   ├── sample1_R2.fastq.gz
│   ├── sample2_R1.fastq.gz
│   └── sample2_R2.fastq.gz
└── README.md         # This file
```

**Note:** The actual data files are ignored by git (see `.gitignore`) to prevent committing large files.

---

## Setup Instructions

### 1. Place Raw FASTQ Files

Copy your raw Illumina NovaSeq paired-end reads into `data/raw/`:

```bash
# Example: Copy files from sequencing run
cp /path/to/sequencing/run/*_R1.fastq.gz data/raw/
cp /path/to/sequencing/run/*_R2.fastq.gz data/raw/
```

**File naming requirements:**
- Paired-end reads must have `_R1` and `_R2` in filenames
- Files should be gzip compressed (`.fastq.gz`)
- Sample names should not contain spaces or special characters

### 2. Update Sample Sheet

Edit `config/samples.tsv` to match your data files:

```tsv
sample	r1	r2
sample1	data/raw/sample1_R1.fastq.gz	data/raw/sample1_R2.fastq.gz
sample2	data/raw/sample2_R1.fastq.gz	data/raw/sample2_R2.fastq.gz
```

Or update `config/config.yaml` directly:

```yaml
samples:
  sample1:
    r1: "data/raw/sample1_R1.fastq.gz"
    r2: "data/raw/sample1_R2.fastq.gz"
  sample2:
    r1: "data/raw/sample2_R1.fastq.gz"
    r2: "data/raw/sample2_R2.fastq.gz"
```

---

## Expected File Specifications

### Sequencing Platform
- **Platform:** Illumina NovaSeq (2-channel chemistry)
- **Read type:** Paired-end
- **Read length:** Typically 150bp

### VLP Sample Requirements
- Samples should be from VLP (Virus-Like Particle) enriched preparations
- RdAB (Random displacement Amplification) protocol
- Library prep: Mechanical shearing or tagmentation-based

### File Format
- **Format:** FASTQ (gzipped)
- **Quality scores:** Phred+33 encoding
- **File size:** Varies (typically 1-5 GB per file)

---

## Troubleshooting

### File Not Found Errors

**Problem:** Pipeline can't find input files

**Solution:**
```bash
# Check files exist
ls -lh data/raw/

# Verify filenames match config
cat config/samples.tsv

# Check file permissions
chmod 644 data/raw/*.fastq.gz
```

### Large File Warnings

**Problem:** Files are very large (>10 GB each)

**Note:** This is normal for high-coverage sequencing. Ensure you have sufficient disk space:
- Raw data: ~10-50 GB
- Pipeline output: ~20-100 GB
- Total required: ~100-200 GB

### Empty Files

**Problem:** Downloaded/transferred files are empty or corrupt

**Solution:**
```bash
# Check file sizes
ls -lh data/raw/

# Verify gzip integrity
gunzip -t data/raw/*.fastq.gz

# Check read counts
zcat data/raw/sample1_R1.fastq.gz | wc -l
# Should be: (number of reads × 4)
```

---

## Test Data

For testing the pipeline without full datasets, you can create subsampled files:

```bash
# Subsample first 100,000 reads for testing
zcat data/raw/sample1_R1.fastq.gz | head -n 400000 | gzip > data/raw/test_R1.fastq.gz
zcat data/raw/sample1_R2.fastq.gz | head -n 400000 | gzip > data/raw/test_R2.fastq.gz
```

Then update config to use test files and run a quick test:

```bash
snakemake --use-conda --cores 4 -n  # dry run
snakemake --use-conda --cores 4      # execute
```

---

## Data Organization Best Practices

### File Naming Convention

Use clear, consistent naming:

```
Good:
- VLP001_S1_R1.fastq.gz
- VLP001_S1_R2.fastq.gz
- Feces_TP1_R1.fastq.gz

Bad:
- file 1.fastq.gz  (spaces)
- sample-final-v2-FINAL.fastq.gz  (unclear)
- R1.fastq.gz  (no sample ID)
```

### Directory Organization

For multiple experiments or time points:

```
data/
├── raw/
│   ├── experiment1/
│   │   ├── sample1_R1.fastq.gz
│   │   └── sample1_R2.fastq.gz
│   └── experiment2/
│       ├── sample2_R1.fastq.gz
│       └── sample2_R2.fastq.gz
```

Then update config paths accordingly.

### Metadata Tracking

Keep a separate metadata file with sample information:

```bash
# data/metadata.tsv
sample_id	date	subject	timepoint	notes
VLP001	2025-01-15	P001	baseline	Good quality
VLP002	2025-01-22	P001	week1	Repeat needed
```

---

## Storage Recommendations

### Local Development
- Keep only test/small datasets locally
- Use external drives for full datasets
- Clean up intermediate files regularly

### HPC/Cluster
- Use scratch/tmp directories for intermediate files
- Archive final results to permanent storage
- Compress old datasets when not in use

### Backups
- Keep raw data backed up in at least 2 locations
- Raw FASTQ files are irreplaceable
- Pipeline outputs can be regenerated

---

## Questions?

- Check pipeline README: `../README.md`
- Review configuration: `../config/config.yaml`
- Open an issue: https://github.com/shandley/lab-virome-QC/issues
- Contact: scott.handley@wustl.edu
