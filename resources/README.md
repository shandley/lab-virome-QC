# Reference Database Resources

This directory contains reference databases required for the Lab Virome QC pipeline.

## Required Files

### 1. PhiX174 Reference Genome
**File:** `phix174.fasta`
**Purpose:** Remove Illumina PhiX control sequences
**Size:** ~5 KB
**Source:** NCBI RefSeq

PhiX174 is a standard Illumina sequencing control that must be removed from virome data.

### 2. Host Genome Reference
**File:** `host_genome.fasta`
**Purpose:** Remove host contamination sequences
**Size:** Varies (human: ~3 GB)
**Source:** Ensembl, NCBI, or UCSC

For VLP-enriched samples, high host contamination (>10%) indicates VLP preparation failure.

**Common host organisms:**
- Human (GRCh38)
- Mouse (GRCm39)
- Other mammalian hosts

### 3. rRNA Database
**File:** `silva_rrna.fasta`
**Purpose:** Remove ribosomal RNA contamination
**Size:** ~100-500 MB
**Source:** SILVA database

Critical for RT-based protocols where all RNA gets converted to cDNA.

---

## Quick Setup

### Automated Setup (Recommended)

Run the setup script from the repository root:

```bash
# Download all required references
bash scripts/setup_references.sh

# For specific host organism:
bash scripts/setup_references.sh --host human
bash scripts/setup_references.sh --host mouse
```

### Manual Setup

If you need to customize or troubleshoot, use individual download scripts:

```bash
# 1. Download PhiX174 reference
bash resources/download_phix.sh

# 2. Download host genome (human example)
bash resources/download_host.sh human

# 3. Download SILVA rRNA database
bash resources/download_silva.sh
```

---

## Detailed Setup Instructions

### PhiX174 Reference

**Option A: From NCBI (Recommended)**
```bash
wget -O resources/phix174.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=fasta&retmode=text"
```

**Option B: From BBTools**
PhiX reference is also included with BBTools installation.

### Host Genome Reference

**For Human (GRCh38):**
```bash
# Download from Ensembl
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa resources/host_genome.fasta

# Create minimap2 index (speeds up mapping)
minimap2 -d resources/host_genome.mmi resources/host_genome.fasta
```

**For Mouse (GRCm39):**
```bash
wget ftp://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
mv Mus_musculus.GRCm39.dna.primary_assembly.fa resources/host_genome.fasta
```

**For Other Organisms:**
- Visit Ensembl: https://www.ensembl.org/
- Or NCBI: https://www.ncbi.nlm.nih.gov/genome/

### SILVA rRNA Database

**Option A: SILVA SSU/LSU (Recommended for viromes)**
```bash
# Download SILVA SSU+LSU database (rRNA only)
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz

# Combine SSU and LSU
gunzip SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
gunzip SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz
cat SILVA_138.1_SSURef_NR99_tax_silva.fasta SILVA_138.1_LSURef_NR99_tax_silva.fasta > resources/silva_rrna.fasta

# Clean up
rm SILVA_138.1_SSURef_NR99_tax_silva.fasta SILVA_138.1_LSURef_NR99_tax_silva.fasta
```

**Option B: BBTools rRNA database**
BBTools includes a pre-built rRNA database that can also be used.

---

## Verification

After downloading, verify all files are present:

```bash
ls -lh resources/
```

Expected output:
```
phix174.fasta           # ~5 KB
host_genome.fasta       # ~3 GB (human)
silva_rrna.fasta        # ~100-500 MB
```

Test that FASTA files are valid:
```bash
# Count sequences in each file
grep -c "^>" resources/phix174.fasta        # Should be 1
grep -c "^>" resources/host_genome.fasta    # Varies by organism
grep -c "^>" resources/silva_rrna.fasta     # Should be thousands
```

---

## Storage Requirements

| File | Organism | Size (compressed) | Size (uncompressed) |
|------|----------|-------------------|---------------------|
| PhiX174 | N/A | ~5 KB | ~5 KB |
| Host genome | Human (GRCh38) | ~900 MB | ~3.1 GB |
| Host genome | Mouse (GRCm39) | ~800 MB | ~2.7 GB |
| SILVA rRNA | N/A | ~50 MB | ~350 MB |

**Total:** ~4-5 GB for human samples

---

## Updating References

### When to Update

- **PhiX174:** Rarely changes (stable reference)
- **Host genome:** Update yearly or when new assembly released
- **SILVA:** Update with each major release (every 1-2 years)

### Version Tracking

Document which versions you're using:

```bash
# Create version file
echo "PhiX174: NC_001422.1" > resources/VERSIONS.txt
echo "Host: GRCh38 (Ensembl release 110)" >> resources/VERSIONS.txt
echo "SILVA: v138.1" >> resources/VERSIONS.txt
echo "Downloaded: $(date)" >> resources/VERSIONS.txt
```

---

## Troubleshooting

### Download Fails

**Issue:** wget/curl timeout or connection issues

**Solution:**
```bash
# Try with curl instead
curl -o resources/phix174.fasta "https://..."

# Or use aria2c for resumable downloads
aria2c -x 4 -s 4 [URL]
```

### File Corrupted

**Issue:** Download interrupted or corrupted file

**Solution:**
```bash
# Verify file integrity if md5sum provided
md5sum resources/host_genome.fasta

# Re-download if needed
rm resources/host_genome.fasta
bash resources/download_host.sh human
```

### Disk Space Issues

**Issue:** Not enough space for large genome files

**Solution:**
```bash
# Keep files compressed when possible
gzip resources/host_genome.fasta
# Update config to use .fasta.gz

# Or use smaller reference (masked genome)
# Or use minimap2 index only (.mmi files)
```

### Permission Denied

**Issue:** Cannot write to resources/ directory

**Solution:**
```bash
# Check permissions
ls -ld resources/

# Fix permissions if needed
chmod 755 resources/
```

---

## Custom References

### Using Your Own Host Genome

If working with a non-standard organism:

```bash
# 1. Place your genome FASTA in resources/
cp /path/to/your_genome.fasta resources/host_genome.fasta

# 2. Update config/config.yaml if needed
# references:
#   host_genome: "resources/host_genome.fasta"

# 3. Index for minimap2 (optional but recommended)
minimap2 -d resources/host_genome.mmi resources/host_genome.fasta
```

### Using Alternative rRNA Databases

```bash
# RDP Ribosomal Database
# Download from: https://rdp.cme.msu.edu/

# Or use your own curated rRNA sequences
# Just ensure FASTA format
```

---

## Citation

If using these databases in publications, please cite:

**SILVA:**
- Quast C et al. (2013) The SILVA ribosomal RNA gene database project. Nucleic Acids Research 41:D590-D596

**Ensembl:**
- Cunningham F et al. (2022) Ensembl 2022. Nucleic Acids Research 50:D988-D995

**RefSeq:**
- O'Leary NA et al. (2016) Reference sequence (RefSeq) database at NCBI. Nucleic Acids Research 44:D733-D745

---

## Questions?

- Open an issue: https://github.com/shandley/lab-virome-QC/issues
- Email: scott.handley@wustl.edu
- Lab Slack: #virome-qc
