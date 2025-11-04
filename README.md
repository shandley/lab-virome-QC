# Lab Virome QC Pipeline

A comprehensive quality control pipeline for VLP-enriched virome sequencing data generated from RdAB (Random displacement Amplification) protocol and Illumina NovaSeq sequencing.

---

## Overview

This Snakemake pipeline provides robust QC specifically designed for:

- **VLP (Virus-Like Particle) enriched samples**
- **RdAB amplification** (RT + Random priming + PCR)
- **Illumina NovaSeq sequencing** (2-channel chemistry)
- **Mechanical shearing or tagmentation-based** library prep

### Key Features

✅ **NovaSeq-specific QC** - PolyG tail removal (critical for 2-channel chemistry)
✅ **VLP enrichment assessment** - ViromeQC enrichment scoring
✅ **Comprehensive contamination removal** - PhiX, host, rRNA
✅ **Optical duplicate removal** - Illumina patterned flow cell artifacts
✅ **Automated QC flagging** - Pass/fail criteria for each sample
✅ **Rich reporting** - MultiQC dashboard with all metrics

---

## Pipeline Workflow

```
Raw Reads (NovaSeq FASTQ)
    ↓
[1] FastQC (raw reads)
    ↓
[2] Clumpify (remove optical duplicates)
    ↓
[3] fastp (CRITICAL: polyG removal + adapter trimming + QC)
    ↓
[4] FastQC (trimmed reads)
    ↓
[5] BBDuk (PhiX removal)
    ↓
[6] minimap2 (host depletion) ← QC metric for VLP success
    ↓
[7] BBDuk (rRNA removal)
    ↓
[8] ViromeQC (enrichment assessment) ← PRIMARY QC metric
    ↓
[9] FastQC (final clean reads)
    ↓
[10] MultiQC (aggregate all reports)
    ↓
Clean reads + QC reports + Sample flags
```

---

## Installation

### Requirements

- [Conda/Mamba](https://github.com/conda-forge/miniforge) (for environment management)
- [Snakemake](https://snakemake.readthedocs.io/) ≥7.0

### Quick Start

```bash
# Clone repository
git clone https://github.com/YOUR_USERNAME/lab-virome-QC.git
cd lab-virome-QC

# Install Snakemake (if not already installed)
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake

# Test installation
snakemake --version
```

---

## Setup

### 1. Prepare Reference Databases

**Quick Setup (Recommended):**

Download all required reference databases with one command:

```bash
bash scripts/setup_references.sh human
```

This will download:
- PhiX174 reference (~5 KB)
- Human genome GRCh38 (~900 MB compressed, ~3 GB uncompressed)
- SILVA rRNA database (~150 MB compressed, ~350 MB uncompressed)

**Time:** 15-30 minutes depending on internet speed
**Disk space:** ~4-5 GB total

**For other organisms:**
```bash
bash scripts/setup_references.sh mouse    # For mouse samples
bash scripts/setup_references.sh custom   # For custom genome (will prompt for URL)
```

**Manual/Individual Downloads:**

If you prefer to download databases individually:

```bash
bash resources/download_phix.sh           # PhiX174 reference
bash resources/download_host.sh human     # Host genome
bash resources/download_silva.sh          # SILVA rRNA database
```

**For detailed documentation, troubleshooting, and advanced options:**
See `resources/README.md`

### 2. Configure Sample Information

Edit `config/config.yaml` or create a sample sheet `config/samples.tsv`:

**Option A: Edit config.yaml directly**
```yaml
samples:
  sample1:
    r1: "data/raw/sample1_R1.fastq.gz"
    r2: "data/raw/sample1_R2.fastq.gz"
  sample2:
    r1: "data/raw/sample2_R1.fastq.gz"
    r2: "data/raw/sample2_R2.fastq.gz"
```

**Option B: Use sample sheet (recommended for many samples)**
```bash
# Edit config/samples.tsv
sample	r1	r2
sample1	data/raw/sample1_R1.fastq.gz	data/raw/sample1_R2.fastq.gz
sample2	data/raw/sample2_R1.fastq.gz	data/raw/sample2_R2.fastq.gz
```

### 3. Adjust QC Thresholds (Optional)

Edit `config/config.yaml` to set custom QC thresholds:

```yaml
qc_thresholds:
  min_enrichment_score: 10      # ViromeQC enrichment score
  max_host_percent: 10          # Maximum % host reads
  max_rrna_percent: 20          # Maximum % rRNA after removal
  min_final_reads: 100000       # Minimum reads after QC
```

---

## Usage

### Run Complete Pipeline

```bash
# Dry run (check what will be executed)
snakemake --use-conda -n

# Run locally with 8 cores
snakemake --use-conda --cores 8

# Run on HPC with SLURM
snakemake --use-conda --profile slurm --jobs 20
```

### Run Specific Steps

```bash
# Just run FastQC on raw reads
snakemake --use-conda --cores 4 results/fastqc/raw/

# Run up to adapter trimming
snakemake --use-conda --cores 8 results/fastp/

# Generate MultiQC report only
snakemake --use-conda --cores 2 results/multiqc/multiqc_report.html
```

### Generate DAG Visualization

```bash
snakemake --dag | dot -Tpng > dag.png
```

---

## Output Structure

```
results/
├── fastqc/                    # FastQC reports (raw, trimmed, final)
├── clumpify/                  # Optical duplicate removal
├── fastp/                     # Adapter trimming + QC
├── phix_removed/              # PhiX-depleted reads
├── host_depleted/             # Host-depleted reads
├── rrna_removed/              # rRNA-depleted reads (clean)
├── viromeqc/                  # ViromeQC enrichment scores
├── clean_reads/               # Symlinks to final clean reads
├── reports/
│   ├── read_counts.tsv        # Read counts at each step
│   └── sample_qc_flags.tsv    # Pass/fail flags per sample
├── multiqc/
│   └── multiqc_report.html    # Comprehensive QC dashboard
└── logs/                      # All log files
```

---

## Interpreting Results

### 1. MultiQC Report

Open `results/multiqc/multiqc_report.html` in a web browser.

**Key sections to check:**
- **FastQC**: Look for polyG in overrepresented sequences (should be gone after fastp)
- **fastp**: Check adapter content, insert size distribution, duplication rates
- **Read counts**: Track reads retained at each step

### 2. Sample QC Flags

Check `results/reports/sample_qc_flags.tsv`:

```
sample    enrichment_score  pass_enrichment  pass_host  pass_rrna  overall_pass  notes
sample1   25.3              PASS             PASS       PASS       PASS          None
sample2   3.2               FAIL             FAIL       PASS       FAIL          Low_enrichment(3.2);High_host(15.3%)
```

**Interpretation:**
- **PASS**: Sample meets all QC criteria
- **FAIL**: Sample fails one or more criteria (check notes)
- **WARNING**: Insufficient data to assess

### 3. ViromeQC Enrichment Score

**Most important QC metric for VLP samples!**

- **Score >>1 (e.g., >10)**: Good VLP enrichment
- **Score ~1**: Poor enrichment (essentially a metagenome)
- **Score <1**: VLP prep failed

Low scores indicate:
- Failed VLP preparation
- Excessive bacterial contamination
- Sample may need to be excluded or re-processed

### 4. Host Contamination

Check `results/host_depleted/*_host_stats.txt`

- **<5% host reads**: Excellent VLP prep
- **5-10% host**: Acceptable
- **>10% host**: VLP prep likely failed

### 5. Read Retention

Check `results/reports/read_counts.tsv`

Typical read retention through QC pipeline:
```
Raw → Clean: 20-50% retention is normal for VLP samples
```

Major losses expected at:
- fastp (adapter trimming, polyG, quality)
- rRNA removal (if RT-based protocol)

---

## Critical NovaSeq + VLP Considerations

### PolyG Tails (NovaSeq Artifact)

**Problem:** NovaSeq 2-channel chemistry calls dark cycles as high-quality G bases
**Solution:** fastp with `--trim_poly_g` (enabled by default in this pipeline)
**Check:** FastQC on raw reads should show GGGG in overrepresented sequences; should be gone after fastp

### VLP Enrichment Assessment

**Problem:** VLP prep can fail, resulting in metagenome-like contamination
**Solution:** ViromeQC quantifies enrichment score
**Action:** Flag samples with enrichment score <10 for review/exclusion

### Host Contamination as QC Metric

**Problem:** High host reads indicate VLP prep failure
**Solution:** Track % host reads; flag >10%
**Note:** Unlike metagenome QC, host contamination is QC failure, not just data cleanup

---

## Troubleshooting

### Pipeline fails at ViromeQC

**Error:** `viromeQC.py: command not found`

**Solution:**
```bash
# Install ViromeQC manually
conda activate viromeqc
pip install viromeqc
```

### High memory usage

**Problem:** BBTools rules consuming too much RAM

**Solution:** Reduce memory allocation in Snakefile:
```python
resources:
    mem_mb = 8000  # Reduce from 16000
```

### Pipeline is slow

**Problem:** Single-threaded execution

**Solution:**
```bash
# Use more cores
snakemake --use-conda --cores 16

# Or run on cluster
snakemake --use-conda --cluster "sbatch" --jobs 50
```

### Adapter sequences not removed

**Problem:** Custom primers not detected by fastp

**Solution:** Add custom adapter sequences to fastp rule:
```bash
--adapter_sequence CUSTOM_PRIMER_B_SEQ \
--adapter_sequence_r2 CUSTOM_PRIMER_B_R2_SEQ
```

---

## Citation

If you use this pipeline in your research, please cite:

- **Snakemake:** Mölder et al. (2021) F1000Research
- **fastp:** Chen et al. (2018) Bioinformatics
- **ViromeQC:** Zolfo et al. (2019) Nature Biotechnology
- **BBTools:** Bushnell (2014) BBMap
- **MultiQC:** Ewels et al. (2016) Bioinformatics

---

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contributing

We welcome contributions from lab members and the community!

### For Lab Members

- **New to the project?** Start with [LEARNING_RESOURCES.md](docs/LEARNING_RESOURCES.md)
- **Ready to contribute?** Read [CONTRIBUTING.md](CONTRIBUTING.md)
- **Need help with GitHub?** See [GITHUB_SETUP.md](docs/GITHUB_SETUP.md)
- **Questions?** Check our [Code of Conduct](CODE_OF_CONDUCT.md)

### Quick Start for Contributors

```bash
# Fork and clone
git clone https://github.com/YOUR_USERNAME/lab-virome-QC.git
cd lab-virome-QC

# Create a branch
git checkout -b feature/my-feature

# Make changes, commit, push
git add .
git commit -m "feat: description of changes"
git push origin feature/my-feature

# Open a Pull Request on GitHub
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed guidelines.

---

## Contact

For questions or issues:
- Open an [Issue](https://github.com/shandley/lab-virome-QC/issues)
- Email: scott.handley@wustl.edu
- Lab Slack: #virome-qc

---

## Acknowledgments

Developed for VLP-enriched virome analysis workflows.

Special considerations for:
- NovaSeq 2-channel chemistry artifacts
- RdAB amplification biases
- VLP enrichment quality assessment
- Mechanical shearing library preparation

---

## References

1. Chen et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 34(17):i884-i890.
2. Zolfo et al. (2019). Detecting contamination in viromes using ViromeQC. *Nature Biotechnology* 37:1408-1412.
3. Bushnell (2014). BBMap: A Fast, Accurate, Splice-Aware Aligner. *Lawrence Berkeley National Lab*
4. Ewels et al. (2016). MultiQC: summarize analysis results for multiple tools and samples. *Bioinformatics* 32(19):3047-3048.
5. 2024 benchmarking studies on virome QC best practices
