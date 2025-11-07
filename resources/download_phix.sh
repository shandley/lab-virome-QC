#!/bin/bash
# Download PhiX174 reference genome for Illumina control sequence removal

set -e  # Exit on error
set -u  # Exit on undefined variable

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}======================================${NC}"
echo -e "${GREEN}PhiX174 Reference Download${NC}"
echo -e "${GREEN}======================================${NC}"

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUTPUT_FILE="${SCRIPT_DIR}/phix174.fasta"

# Check if file already exists
if [ -f "$OUTPUT_FILE" ]; then
    echo -e "${YELLOW}Warning: phix174.fasta already exists${NC}"
    read -p "Do you want to overwrite it? (y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo -e "${GREEN}Keeping existing file. Exiting.${NC}"
        exit 0
    fi
    rm "$OUTPUT_FILE"
fi

echo -e "\n${GREEN}Downloading PhiX174 reference from NCBI...${NC}"

# Download from NCBI using Entrez utilities
# NC_001422.1 is the PhiX174 RefSeq accession
if command -v wget &> /dev/null; then
    wget -O "$OUTPUT_FILE" \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=fasta&retmode=text" \
        2>&1 | grep -v "^--" || true
elif command -v curl &> /dev/null; then
    curl -o "$OUTPUT_FILE" \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=fasta&retmode=text"
else
    echo -e "${RED}Error: Neither wget nor curl is installed${NC}"
    echo "Please install wget or curl and try again"
    exit 1
fi

# Verify download
if [ ! -f "$OUTPUT_FILE" ]; then
    echo -e "${RED}Error: Download failed${NC}"
    exit 1
fi

# Check file is not empty
if [ ! -s "$OUTPUT_FILE" ]; then
    echo -e "${RED}Error: Downloaded file is empty${NC}"
    rm "$OUTPUT_FILE"
    exit 1
fi

# Verify FASTA format
if ! grep -q "^>" "$OUTPUT_FILE"; then
    echo -e "${RED}Error: File does not appear to be in FASTA format${NC}"
    rm "$OUTPUT_FILE"
    exit 1
fi

# Get file info
FILE_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
SEQ_COUNT=$(grep -c "^>" "$OUTPUT_FILE")

echo -e "\n${GREEN}======================================${NC}"
echo -e "${GREEN}Download Complete!${NC}"
echo -e "${GREEN}======================================${NC}"
echo -e "File: ${OUTPUT_FILE}"
echo -e "Size: ${FILE_SIZE}"
echo -e "Sequences: ${SEQ_COUNT}"
echo -e ""

# Display first few lines
echo -e "${GREEN}First lines of file:${NC}"
head -n 5 "$OUTPUT_FILE"
echo ""

echo -e "${GREEN}PhiX174 reference is ready to use!${NC}"
