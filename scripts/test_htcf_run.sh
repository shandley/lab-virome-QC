#!/bin/bash
# =============================================================================
# HTCF Test Run Script for Lab Virome QC Pipeline
# =============================================================================
#
# This script performs a test run of the pipeline on HTCF
# Use this to verify the pipeline works before running on full dataset
#
# Usage:
#   bash scripts/test_htcf_run.sh
#
# Prerequisites:
#   1. Connected to HTCF: ssh <username>@login.htcf.wustl.edu
#   2. Snakemake environment activated
#   3. Reference databases downloaded
#   4. Test data available in data/raw/
#
# =============================================================================

set -euo pipefail  # Exit on error, undefined variable, or pipe failure

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# =============================================================================
# Functions
# =============================================================================

print_header() {
    echo -e "\n${BLUE}======================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}======================================${NC}\n"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

print_info() {
    echo -e "${BLUE}ℹ $1${NC}"
}

# =============================================================================
# Pre-flight Checks
# =============================================================================

print_header "HTCF Pipeline Test - Pre-flight Checks"

# Check if we're on HTCF
if [[ ! $(hostname) =~ htcf ]]; then
    print_warning "Not running on HTCF. This script is designed for HTCF."
    read -p "Continue anyway? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# Check if Snakemake is available
if command -v snakemake &> /dev/null; then
    SNAKEMAKE_VERSION=$(snakemake --version)
    print_success "Snakemake found (version $SNAKEMAKE_VERSION)"
else
    print_error "Snakemake not found!"
    echo "Please activate the Snakemake environment:"
    echo "  source /ref/sahlab/software/miniforge3/bin/activate"
    echo "  conda activate snakemake_tutorial"
    exit 1
fi

# Check if SLURM commands are available
if command -v sbatch &> /dev/null; then
    print_success "SLURM commands available"
else
    print_error "SLURM commands not found!"
    echo "Make sure you're on the HTCF login node, not a compute node."
    exit 1
fi

# Check if config file exists
if [[ -f config/config.yaml ]]; then
    print_success "Configuration file found"
else
    print_error "config/config.yaml not found!"
    exit 1
fi

# Check if SLURM profile exists
if [[ -f profiles/htcf/config.yaml ]]; then
    print_success "HTCF SLURM profile found"
else
    print_error "profiles/htcf/config.yaml not found!"
    echo "The SLURM profile is required to run on HTCF."
    exit 1
fi

# Check for reference databases
print_info "Checking reference databases..."
MISSING_REFS=0

if [[ ! -f resources/phix174.fasta ]]; then
    print_warning "PhiX reference not found"
    MISSING_REFS=1
fi

if [[ ! -f resources/host_genome.fasta ]]; then
    print_warning "Host genome not found"
    MISSING_REFS=1
fi

if [[ ! -f resources/silva_rrna.fasta ]]; then
    print_warning "SILVA rRNA database not found"
    MISSING_REFS=1
fi

if [[ $MISSING_REFS -eq 1 ]]; then
    print_warning "Some reference databases are missing"
    echo "Download them with:"
    echo "  bash scripts/setup_references.sh human"
    read -p "Continue anyway? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
else
    print_success "All reference databases found"
fi

# Check for test data
print_info "Checking for test data..."
FASTQ_COUNT=$(find data/raw -name "*.fastq.gz" 2>/dev/null | wc -l)

if [[ $FASTQ_COUNT -eq 0 ]]; then
    print_error "No FASTQ files found in data/raw/"
    echo "Please add test data to data/raw/ before running this script"
    exit 1
else
    print_success "Found $FASTQ_COUNT FASTQ file(s) in data/raw/"
fi

# Create SLURM log directory
mkdir -p logs/slurm
print_success "SLURM log directory created"

# =============================================================================
# Dry Run
# =============================================================================

print_header "Step 1: Dry Run (Checking Pipeline Logic)"

print_info "Running Snakemake dry run to check for errors..."

if snakemake --profile profiles/htcf --dry-run; then
    print_success "Dry run completed successfully!"
    echo ""
    snakemake --profile profiles/htcf --dry-run --quiet 2>&1 | grep -E "Job stats:|rule" | tail -20
else
    print_error "Dry run failed!"
    echo "Please check the error messages above and fix configuration issues."
    exit 1
fi

# =============================================================================
# Ask User to Proceed
# =============================================================================

echo ""
read -p "$(echo -e ${YELLOW}Dry run successful. Proceed with actual test run? (y/n)${NC}) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    print_info "Test run cancelled by user"
    exit 0
fi

# =============================================================================
# Test Run
# =============================================================================

print_header "Step 2: Running Pipeline on HTCF"

print_info "Submitting pipeline jobs to SLURM..."
print_info "This may take several hours depending on data size"
print_info ""
print_info "Monitor progress with:"
print_info "  - squeue -u \$USER                  (check job queue)"
print_info "  - tail -f .snakemake/log/*.log     (view Snakemake log)"
print_info "  - tail -f logs/slurm/*.out         (view individual job logs)"
echo ""

# Run with limited concurrent jobs for testing
snakemake --profile profiles/htcf --jobs 20

# =============================================================================
# Check Results
# =============================================================================

print_header "Step 3: Checking Results"

# Check if MultiQC report was generated
if [[ -f results/multiqc/multiqc_report.html ]]; then
    print_success "MultiQC report generated"
else
    print_error "MultiQC report not found!"
fi

# Check if QC flags were generated
if [[ -f results/reports/sample_qc_flags.tsv ]]; then
    print_success "QC flags generated"
    echo ""
    echo "Sample QC Summary:"
    echo "=================="
    cat results/reports/sample_qc_flags.tsv
else
    print_error "QC flags file not found!"
fi

# Check if read counts were generated
if [[ -f results/reports/read_counts.tsv ]]; then
    print_success "Read counts generated"
else
    print_error "Read counts file not found!"
fi

# Check for clean reads
CLEAN_READS=$(find results/clean_reads -name "*.fastq.gz" 2>/dev/null | wc -l)
if [[ $CLEAN_READS -gt 0 ]]; then
    print_success "Found $CLEAN_READS clean read file(s)"
else
    print_warning "No clean reads found"
fi

# =============================================================================
# Summary
# =============================================================================

print_header "Test Run Complete!"

echo "Next steps:"
echo "==========="
echo "1. View the MultiQC report:"
echo "   - Download: scp <user>@login.htcf.wustl.edu:$(pwd)/results/multiqc/multiqc_report.html ."
echo "   - Open in browser: open multiqc_report.html"
echo ""
echo "2. Check sample QC flags:"
echo "   - cat results/reports/sample_qc_flags.tsv"
echo ""
echo "3. Review logs for any warnings:"
echo "   - ls -lh results/logs/"
echo ""
echo "4. If test successful, run on full dataset:"
echo "   - snakemake --profile profiles/htcf --jobs 50"
echo ""

print_success "Test completed successfully!"
