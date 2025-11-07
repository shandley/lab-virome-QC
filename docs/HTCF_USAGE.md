# Running Lab Virome QC on HTCF (WashU)

This guide explains how to run the Lab Virome QC pipeline on Washington University's High Throughput Computing Facility (HTCF).

> **Using a different cluster?** This guide is specific to WashU HTCF, but the configuration serves as a working example. See [Adapting for Your Own Cluster](#adapting-for-your-own-cluster) to learn how to modify it for your institution.

---

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Initial Setup](#initial-setup)
3. [Prepare Your Data](#prepare-your-data)
4. [Configure the Pipeline](#configure-the-pipeline)
5. [Run the Pipeline](#run-the-pipeline)
6. [Monitor Jobs](#monitor-jobs)
7. [Retrieve Results](#retrieve-results)
8. [Troubleshooting](#troubleshooting)
9. [Adapting for Your Own Cluster](#adapting-for-your-own-cluster)

---

## Prerequisites

- Active HTCF account
- SSH access to HTCF: `ssh <your_username>@login.htcf.wustl.edu`
- Sequencing data (FASTQ files) uploaded to HTCF
- Basic familiarity with Linux command line

---

## Initial Setup

### 1. Connect to HTCF

```bash
ssh leranwang@login.htcf.wustl.edu
```

### 2. Clone the Repository

```bash
# Navigate to your working directory
cd /scratch/sahlab/<your_username>  # or your preferred location

# Clone the repository
git clone https://github.com/shandley/lab-virome-QC.git
cd lab-virome-QC
```

### 3. Activate Snakemake Environment

```bash
# Activate the lab's Snakemake environment
source /ref/sahlab/software/miniforge3/bin/activate
conda activate snakemake_tutorial

# Verify Snakemake is available
snakemake --version
```

### 4. Download Reference Databases

**Important:** Only do this once! References can be shared across projects.

```bash
# Option 1: Download to a shared location (recommended)
# Check with lab if references already exist in /ref/sahlab/

# Option 2: Download to your project directory
bash scripts/setup_references.sh human

# This will download (~4-5 GB total):
# - PhiX174 reference
# - Human genome GRCh38
# - SILVA rRNA database
```

**Time:** 15-30 minutes depending on network speed

---

## Prepare Your Data

### 1. Organize Your FASTQ Files

Create the data directory structure:

```bash
mkdir -p data/raw
```

### 2. Link or Copy Your FASTQ Files

**Option A: Symbolic Links (Recommended - saves space)**

```bash
# Link your existing data
ln -s /path/to/your/data/*_R1*.fastq.gz data/raw/
ln -s /path/to/your/data/*_R2*.fastq.gz data/raw/
```

**Option B: Copy Files**

```bash
# Copy data to project directory
cp /path/to/your/data/*.fastq.gz data/raw/
```

### 3. Verify Files

```bash
ls -lh data/raw/
```

Expected format:
```
sample1_R1.fastq.gz
sample1_R2.fastq.gz
sample2_R1.fastq.gz
sample2_R2.fastq.gz
...
```

---

## Configure the Pipeline

### 1. Edit Sample Configuration

**Option A: Edit `config/samples.tsv` (Recommended for many samples)**

```bash
# Edit the sample sheet
nano config/samples.tsv
```

Format (tab-separated):
```
sample	r1	r2
sample1	data/raw/sample1_R1.fastq.gz	data/raw/sample1_R2.fastq.gz
sample2	data/raw/sample2_R1.fastq.gz	data/raw/sample2_R2.fastq.gz
```

**Option B: Edit `config/config.yaml` directly**

```bash
nano config/config.yaml
```

Update the `samples:` section with your sample names and paths.

### 2. Adjust QC Thresholds (Optional)

In `config/config.yaml`, you can customize:

```yaml
qc_thresholds:
  min_enrichment_score: 10      # ViromeQC enrichment score
  max_host_percent: 10          # Maximum % host reads
  max_rrna_percent: 20          # Maximum % rRNA
  min_final_reads: 100000       # Minimum reads after QC
```

### 3. Create SLURM Log Directory

```bash
mkdir -p logs/slurm
```

---

## Run the Pipeline

### Dry Run (Recommended First Step)

Check what Snakemake will execute without actually running anything:

```bash
snakemake --profile profiles/htcf --dry-run
```

This will show:
- All jobs that will be submitted
- Number of jobs per rule
- Expected outputs

### Test Run (Small Dataset)

For initial testing, run on just one sample:

```bash
snakemake --profile profiles/htcf \
    --jobs 10 \
    results/multiqc/multiqc_report.html \
    --config samples="{'sample1': {'r1': 'data/raw/sample1_R1.fastq.gz', 'r2': 'data/raw/sample1_R2.fastq.gz'}}"
```

### Full Pipeline Run

Once the test succeeds, run the full pipeline:

```bash
snakemake --profile profiles/htcf --jobs 50
```

**Parameters:**
- `--profile profiles/htcf`: Use HTCF SLURM configuration
- `--jobs 50`: Submit up to 50 jobs simultaneously

### Alternative: Submit as a Screen Session

For long-running pipelines, use `screen` to keep the process running if you disconnect:

```bash
# Start a screen session
screen -S virome_qc

# Activate environment
source /ref/sahlab/software/miniforge3/bin/activate
conda activate snakemake_tutorial

# Run pipeline
snakemake --profile profiles/htcf --jobs 50

# Detach from screen: Press Ctrl+A, then D

# Reattach later
screen -r virome_qc
```

---

## Monitor Jobs

### Check SLURM Queue

```bash
# View your running jobs
squeue -u $USER

# View all jobs with details
squeue -u $USER -o "%.18i %.9P %.30j %.8u %.2t %.10M %.6D %R"
```

### Check Snakemake Logs

```bash
# View main Snakemake log
tail -f .snakemake/log/*.log

# View SLURM job logs
ls logs/slurm/
tail -f logs/slurm/fastp_sample1.*.out
```

### Check Specific Rule Logs

```bash
# View logs for a specific rule
tail -f results/logs/fastp/sample1.log
```

---

## Retrieve Results

### Check Output Structure

```bash
ls -lh results/
```

Key outputs:
```
results/
├── multiqc/multiqc_report.html         # Main QC dashboard
├── reports/sample_qc_flags.tsv         # Pass/fail flags
├── reports/read_counts.tsv             # Read counts per step
├── viromeqc/                           # ViromeQC scores
├── clean_reads/                        # Final clean FASTQ files
└── logs/                               # All processing logs
```

### Download Results to Local Computer

**From your local computer:**

```bash
# Download MultiQC report
scp leranwang@login.htcf.wustl.edu:/path/to/lab-virome-QC/results/multiqc/multiqc_report.html .

# Download QC flags
scp leranwang@login.htcf.wustl.edu:/path/to/lab-virome-QC/results/reports/sample_qc_flags.tsv .

# Download all results (if small enough)
scp -r leranwang@login.htcf.wustl.edu:/path/to/lab-virome-QC/results/ .
```

---

## Troubleshooting

### Pipeline Fails to Start

**Error:** `Command 'sbatch' not found`

**Solution:** Make sure you're on the login node, not a compute node. Exit any interactive sessions and retry from the login node.

---

### Jobs Stay in PENDING State

**Cause:** SLURM queue is busy or resource requests are too high

**Solution:**
1. Check queue status: `squeue`
2. Reduce concurrent jobs: `snakemake --profile profiles/htcf --jobs 10`
3. Check if you have enough allocation hours

---

### Out of Memory Errors

**Error:** `OUT_OF_MEMORY` in SLURM logs

**Solution:** Edit `profiles/htcf/config.yaml` to increase memory:

```yaml
default-resources:
  - mem_mb=16000  # Increase from 8000
```

Or edit specific rules in `workflow/Snakefile`:

```python
resources:
    mem_mb = 32000  # For memory-intensive rules
```

---

### Conda Environment Creation Fails

**Error:** `CondaEnvException` or package conflicts

**Solution:**
```bash
# Clean conda cache
conda clean --all

# Retry with mamba (faster, better dependency resolution)
snakemake --profile profiles/htcf --use-conda --conda-frontend mamba
```

---

### Reference Databases Missing

**Error:** `FileNotFoundError: resources/host_genome.fasta`

**Solution:**
```bash
# Download missing references
bash scripts/setup_references.sh human

# Or use existing shared references (ask lab)
ln -s /ref/sahlab/references/GRCh38.fa resources/host_genome.fasta
```

---

### Pipeline Stops Unexpectedly

**Cause:** Network disconnection or login session timeout

**Solution:** Use `screen` or `tmux` to maintain persistent sessions:

```bash
# Start screen session
screen -S virome_qc

# Run pipeline inside screen
snakemake --profile profiles/htcf --jobs 50

# Detach: Ctrl+A, then D
# Reattach: screen -r virome_qc
```

---

### Rerunning Failed Jobs

If some jobs fail, Snakemake will automatically rerun only the failed jobs:

```bash
# Rerun with same command
snakemake --profile profiles/htcf --jobs 50
```

To force rerun specific rules:

```bash
# Rerun all fastp jobs
snakemake --profile profiles/htcf --jobs 50 --forcerun fastp_trim
```

---

## Resource Allocation Guidelines

### For Typical Virome Samples (1-5M reads):

| Rule | Recommended Memory | Threads | Time |
|------|-------------------|---------|------|
| fastqc | 4 GB | 2 | 30 min |
| clumpify | 16 GB | 8 | 1 hour |
| fastp | 8 GB | 8 | 1 hour |
| bbduk | 8 GB | 4 | 30 min |
| minimap2 | 16 GB | 8 | 2 hours |
| viromeqc | 8 GB | 4 | 1 hour |
| multiqc | 4 GB | 1 | 10 min |

### For Large Datasets (>10M reads):

Increase memory by 2x and time by 1.5x.

---

## Adapting for Your Own Cluster

**Using a different HPC cluster?** The HTCF profile serves as a working example that you can adapt for your institution's SLURM cluster.

### What to Change

#### 1. Copy the Profile Directory

```bash
# Copy HTCF profile as a starting point
cp -r profiles/htcf profiles/my-cluster
```

#### 2. Edit `profiles/my-cluster/config.yaml`

**Partition Names:**
```yaml
# HTCF uses: general, high-mem, interactive
# Change to your cluster's partitions
default-resources:
  - partition=compute  # Replace with your default partition name
```

Check your available partitions:
```bash
sinfo  # Shows partition names and availability
```

**Memory and Time Defaults:**
```yaml
default-resources:
  - partition=compute      # Your cluster's partition
  - mem_mb=8000           # Adjust based on your cluster's nodes
  - time_min=120          # Adjust based on your cluster's limits
  - tmpdir="/scratch"     # Your cluster's scratch directory
```

**Job Submission Command:**

The `cluster:` line may need adjustment for your scheduler configuration:
```yaml
# If your cluster requires different sbatch parameters:
cluster: "sbatch --partition={resources.partition} --cpus-per-task={threads} --mem={resources.mem_mb}M --time={resources.time_min} --job-name={rule}.{wildcards} --output=logs/slurm/{rule}_{wildcards}.%j.out --error=logs/slurm/{rule}_{wildcards}.%j.err --account=YOUR_ACCOUNT"
```

Common additions:
- `--account=PROJECT_NAME` (if required)
- `--qos=normal` (if using QOS)
- `--constraint=haswell` (if specifying CPU architecture)

#### 3. Test Your Configuration

```bash
# Dry run to check configuration
snakemake --profile profiles/my-cluster --dry-run

# Test with one sample
snakemake --profile profiles/my-cluster --jobs 5 results/fastqc/raw/
```

### Cluster-Specific Considerations

**Check your cluster's documentation for:**

1. **Partition names and limits**
   ```bash
   sinfo -o "%P %l %D %N"
   ```

2. **Memory per node**
   ```bash
   sinfo -o "%P %m %N"
   ```

3. **Account/project requirements**
   ```bash
   sacctmgr show user $USER format=user,account%20
   ```

4. **Time limits per partition**
   ```bash
   scontrol show partition PARTITION_NAME
   ```

### Example Configurations

**Generic University Cluster:**
```yaml
cluster: "sbatch --partition=normal --cpus-per-task={threads} --mem={resources.mem_mb}M --time={resources.time_min} --account=lab_project --job-name={rule}.{wildcards} --output=logs/slurm/{rule}_{wildcards}.%j.out --error=logs/slurm/{rule}_{wildcards}.%j.err"

default-resources:
  - partition=normal
  - mem_mb=8000
  - time_min=240
  - tmpdir="/scratch/$USER"
```

**NIH Biowulf:**
```yaml
cluster: "sbatch --partition=norm --cpus-per-task={threads} --mem={resources.mem_mb}M --time={resources.time_min} --gres=lscratch:10 --job-name={rule}.{wildcards} --output=logs/slurm/{rule}_{wildcards}.%j.out --error=logs/slurm/{rule}_{wildcards}.%j.err"

default-resources:
  - partition=norm
  - mem_mb=8000
  - time_min=240
```

### Need Help?

If you successfully adapt this for your cluster, consider:
- Opening a PR to add your profile as another example
- Sharing your configuration in a GitHub issue to help others
- Adding notes to this documentation

---

## Additional Resources

- **HTCF Documentation:** https://htcf.wustl.edu/docs/
- **Snakemake Documentation:** https://snakemake.readthedocs.io/
- **Lab Slack:** #virome-qc channel
- **Pipeline Issues:** https://github.com/shandley/lab-virome-QC/issues

---

## Getting Help

If you encounter issues:

1. Check SLURM logs: `logs/slurm/`
2. Check rule-specific logs: `results/logs/`
3. Search existing GitHub issues
4. Ask in lab Slack: #virome-qc
5. Open a GitHub issue with:
   - Error message
   - Relevant log files
   - Sample configuration
   - SLURM job ID

---

## Quick Reference Commands

```bash
# Activate environment
source /ref/sahlab/software/miniforge3/bin/activate && conda activate snakemake_tutorial

# Dry run
snakemake --profile profiles/htcf -n

# Run pipeline
snakemake --profile profiles/htcf --jobs 50

# Check jobs
squeue -u $USER

# View logs
tail -f .snakemake/log/*.log

# Cancel all jobs
scancel -u $USER

# Reattach to screen
screen -r virome_qc
```

---

**Last Updated:** November 2025
**Maintainer:** Lab Virome QC Team
