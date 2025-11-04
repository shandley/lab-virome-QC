#!/bin/bash
# Download host genome reference for VLP depletion

set -e  # Exit on error
set -u  # Exit on undefined variable

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${GREEN}======================================${NC}"
echo -e "${GREEN}Host Genome Reference Download${NC}"
echo -e "${GREEN}======================================${NC}"

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUTPUT_FILE="${SCRIPT_DIR}/host_genome.fasta"

# Parse command line arguments
HOST_ORGANISM="${1:-}"

# Display usage if no argument provided
if [ -z "$HOST_ORGANISM" ]; then
    echo -e "${BLUE}Usage: $0 <organism>${NC}"
    echo -e "\nSupported organisms:"
    echo -e "  human    - Homo sapiens (GRCh38)"
    echo -e "  mouse    - Mus musculus (GRCm39)"
    echo -e "  custom   - Provide your own URL"
    echo -e "\nExample:"
    echo -e "  $0 human"
    exit 1
fi

# Check if file already exists
if [ -f "$OUTPUT_FILE" ]; then
    echo -e "${YELLOW}Warning: host_genome.fasta already exists${NC}"
    read -p "Do you want to overwrite it? (y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo -e "${GREEN}Keeping existing file. Exiting.${NC}"
        exit 0
    fi
    rm "$OUTPUT_FILE"
fi

# Set download URL based on organism
case "$HOST_ORGANISM" in
    human)
        echo -e "\n${GREEN}Downloading Human genome (GRCh38) from Ensembl...${NC}"
        echo -e "${YELLOW}Note: This is a large file (~900 MB compressed, ~3 GB uncompressed)${NC}"
        echo -e "${YELLOW}Download may take 10-30 minutes depending on connection speed${NC}\n"

        GENOME_URL="ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        TEMP_FILE="${SCRIPT_DIR}/host_genome.fa.gz"
        ;;

    mouse)
        echo -e "\n${GREEN}Downloading Mouse genome (GRCm39) from Ensembl...${NC}"
        echo -e "${YELLOW}Note: This is a large file (~800 MB compressed, ~2.7 GB uncompressed)${NC}"
        echo -e "${YELLOW}Download may take 10-30 minutes depending on connection speed${NC}\n"

        GENOME_URL="ftp://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
        TEMP_FILE="${SCRIPT_DIR}/host_genome.fa.gz"
        ;;

    custom)
        echo -e "\n${BLUE}Custom genome download${NC}"
        read -p "Enter the download URL: " GENOME_URL
        read -p "Is this a gzipped file? (y/N) " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            TEMP_FILE="${SCRIPT_DIR}/host_genome.fa.gz"
        else
            TEMP_FILE="$OUTPUT_FILE"
        fi
        ;;

    *)
        echo -e "${RED}Error: Unknown organism '$HOST_ORGANISM'${NC}"
        echo -e "Supported: human, mouse, custom"
        exit 1
        ;;
esac

# Download the genome
echo -e "${GREEN}Starting download...${NC}"
if command -v wget &> /dev/null; then
    wget -O "$TEMP_FILE" "$GENOME_URL" || {
        echo -e "${RED}Download failed${NC}"
        rm -f "$TEMP_FILE"
        exit 1
    }
elif command -v curl &> /dev/null; then
    curl -o "$TEMP_FILE" "$GENOME_URL" || {
        echo -e "${RED}Download failed${NC}"
        rm -f "$TEMP_FILE"
        exit 1
    }
else
    echo -e "${RED}Error: Neither wget nor curl is installed${NC}"
    echo "Please install wget or curl and try again"
    exit 1
fi

# Decompress if needed
if [[ "$TEMP_FILE" == *.gz ]]; then
    echo -e "\n${GREEN}Decompressing genome file...${NC}"
    gunzip -c "$TEMP_FILE" > "$OUTPUT_FILE"
    rm "$TEMP_FILE"
fi

# Verify download
if [ ! -f "$OUTPUT_FILE" ]; then
    echo -e "${RED}Error: Final file not created${NC}"
    exit 1
fi

# Check file is not empty
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
echo -e "Sequences (chromosomes/contigs): ${SEQ_COUNT}"
echo -e ""

# Display first few lines
echo -e "${GREEN}First lines of file:${NC}"
head -n 5 "$OUTPUT_FILE"
echo ""

# Optional: Create minimap2 index
echo -e "${YELLOW}Optional: Create minimap2 index for faster mapping?${NC}"
echo -e "This will create a .mmi index file (takes 5-10 minutes)"
read -p "Create index now? (y/N) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    if command -v minimap2 &> /dev/null; then
        echo -e "${GREEN}Creating minimap2 index...${NC}"
        minimap2 -d "${SCRIPT_DIR}/host_genome.mmi" "$OUTPUT_FILE"
        echo -e "${GREEN}Index created: host_genome.mmi${NC}"
    else
        echo -e "${YELLOW}minimap2 not found. Install it to create the index later:${NC}"
        echo -e "conda install -c bioconda minimap2"
        echo -e "minimap2 -d resources/host_genome.mmi resources/host_genome.fasta"
    fi
fi

echo -e "\n${GREEN}Host genome reference is ready to use!${NC}"
