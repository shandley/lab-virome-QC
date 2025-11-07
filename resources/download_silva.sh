#!/bin/bash
# Download SILVA rRNA database for ribosomal RNA removal

set -e  # Exit on error
set -u  # Exit on undefined variable

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${GREEN}======================================${NC}"
echo -e "${GREEN}SILVA rRNA Database Download${NC}"
echo -e "${GREEN}======================================${NC}"

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUTPUT_FILE="${SCRIPT_DIR}/silva_rrna.fasta"

# Check if file already exists
if [ -f "$OUTPUT_FILE" ]; then
    echo -e "${YELLOW}Warning: silva_rrna.fasta already exists${NC}"
    read -p "Do you want to overwrite it? (y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo -e "${GREEN}Keeping existing file. Exiting.${NC}"
        exit 0
    fi
    rm "$OUTPUT_FILE"
fi

# SILVA database URLs (release 138.1)
SILVA_VERSION="138.1"
SSU_URL="https://www.arb-silva.de/fileadmin/silva_databases/release_${SILVA_VERSION//./_}/Exports/SILVA_${SILVA_VERSION}_SSURef_NR99_tax_silva.fasta.gz"
LSU_URL="https://www.arb-silva.de/fileadmin/silva_databases/release_${SILVA_VERSION//./_}/Exports/SILVA_${SILVA_VERSION}_LSURef_NR99_tax_silva.fasta.gz"

echo -e "\n${BLUE}SILVA Database Information:${NC}"
echo -e "Version: ${SILVA_VERSION}"
echo -e "Components: SSU (16S) + LSU (23S) rRNA"
echo -e "Size: ~100-150 MB compressed, ~350 MB uncompressed"
echo -e ""

# Temporary files
SSU_FILE="${SCRIPT_DIR}/silva_ssu.fasta.gz"
LSU_FILE="${SCRIPT_DIR}/silva_lsu.fasta.gz"
SSU_FASTA="${SCRIPT_DIR}/silva_ssu.fasta"
LSU_FASTA="${SCRIPT_DIR}/silva_lsu.fasta"

# Download SSU (16S) rRNA database
echo -e "${GREEN}[1/2] Downloading SSU (16S) rRNA database...${NC}"
if command -v wget &> /dev/null; then
    wget -O "$SSU_FILE" "$SSU_URL" || {
        echo -e "${RED}SSU download failed${NC}"
        rm -f "$SSU_FILE"
        exit 1
    }
elif command -v curl &> /dev/null; then
    curl -L -o "$SSU_FILE" "$SSU_URL" || {
        echo -e "${RED}SSU download failed${NC}"
        rm -f "$SSU_FILE"
        exit 1
    }
else
    echo -e "${RED}Error: Neither wget nor curl is installed${NC}"
    echo "Please install wget or curl and try again"
    exit 1
fi

# Download LSU (23S) rRNA database
echo -e "\n${GREEN}[2/2] Downloading LSU (23S) rRNA database...${NC}"
if command -v wget &> /dev/null; then
    wget -O "$LSU_FILE" "$LSU_URL" || {
        echo -e "${RED}LSU download failed${NC}"
        rm -f "$LSU_FILE" "$SSU_FILE"
        exit 1
    }
elif command -v curl &> /dev/null; then
    curl -L -o "$LSU_FILE" "$LSU_URL" || {
        echo -e "${RED}LSU download failed${NC}"
        rm -f "$LSU_FILE" "$SSU_FILE"
        exit 1
    }
fi

# Decompress files
echo -e "\n${GREEN}Decompressing files...${NC}"
gunzip -c "$SSU_FILE" > "$SSU_FASTA"
gunzip -c "$LSU_FILE" > "$LSU_FASTA"

# Combine SSU and LSU databases
echo -e "${GREEN}Combining SSU and LSU databases...${NC}"
cat "$SSU_FASTA" "$LSU_FASTA" > "$OUTPUT_FILE"

# Clean up temporary files
echo -e "${GREEN}Cleaning up temporary files...${NC}"
rm -f "$SSU_FILE" "$LSU_FILE" "$SSU_FASTA" "$LSU_FASTA"

# Verify final file
if [ ! -f "$OUTPUT_FILE" ]; then
    echo -e "${RED}Error: Final file not created${NC}"
    exit 1
fi

if [ ! -s "$OUTPUT_FILE" ]; then
    echo -e "${RED}Error: File is empty${NC}"
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
echo -e "rRNA Sequences: ${SEQ_COUNT}"
echo -e "Version: SILVA ${SILVA_VERSION}"
echo -e ""

# Display first few lines
echo -e "${GREEN}First lines of file:${NC}"
head -n 3 "$OUTPUT_FILE"
echo ""

# Show statistics
SSU_COUNT=$(grep -c "SSU" "$OUTPUT_FILE" || echo "0")
LSU_COUNT=$(grep -c "LSU" "$OUTPUT_FILE" || echo "0")
echo -e "${BLUE}Database composition:${NC}"
echo -e "  SSU (16S) sequences: ~${SSU_COUNT}"
echo -e "  LSU (23S) sequences: ~${LSU_COUNT}"
echo -e ""

echo -e "${GREEN}SILVA rRNA database is ready to use!${NC}"
echo -e "\n${YELLOW}Citation:${NC}"
echo -e "Quast C et al. (2013) The SILVA ribosomal RNA gene database project."
echo -e "Nucleic Acids Research 41:D590-D596"
