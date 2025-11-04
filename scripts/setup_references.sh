#!/bin/bash
# Master script to set up all reference databases for Lab Virome QC pipeline

set -e  # Exit on error
set -u  # Exit on undefined variable

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Lab Virome QC - Reference Setup${NC}"
echo -e "${GREEN}========================================${NC}"

# Get script directory and project root
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
RESOURCES_DIR="${PROJECT_ROOT}/resources"

# Parse command line arguments
HOST_ORGANISM="${1:-}"
SKIP_PHIX=false
SKIP_HOST=false
SKIP_SILVA=false

# Display usage
usage() {
    echo -e "${BLUE}Usage: $0 [OPTIONS] <host_organism>${NC}"
    echo -e "\nArguments:"
    echo -e "  host_organism    Host organism: human, mouse, or custom"
    echo -e "\nOptions:"
    echo -e "  --skip-phix      Skip PhiX174 download"
    echo -e "  --skip-host      Skip host genome download"
    echo -e "  --skip-silva     Skip SILVA rRNA download"
    echo -e "  -h, --help       Show this help message"
    echo -e "\nExamples:"
    echo -e "  $0 human                    # Download all references (human host)"
    echo -e "  $0 mouse                    # Download all references (mouse host)"
    echo -e "  $0 --skip-host human        # Skip host genome, download others"
    exit 1
}

# Parse options
while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-phix)
            SKIP_PHIX=true
            shift
            ;;
        --skip-host)
            SKIP_HOST=true
            shift
            ;;
        --skip-silva)
            SKIP_SILVA=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            if [ -z "$HOST_ORGANISM" ]; then
                HOST_ORGANISM="$1"
            fi
            shift
            ;;
    esac
done

# Validate host organism if not skipping
if [ "$SKIP_HOST" = false ] && [ -z "$HOST_ORGANISM" ]; then
    echo -e "${RED}Error: Host organism required${NC}\n"
    usage
fi

# Check if resources directory exists
if [ ! -d "$RESOURCES_DIR" ]; then
    echo -e "${RED}Error: resources/ directory not found${NC}"
    echo -e "Please run this script from the project root"
    exit 1
fi

# Summary
echo -e "\n${BLUE}Setup Configuration:${NC}"
echo -e "Project root: ${PROJECT_ROOT}"
echo -e "Resources directory: ${RESOURCES_DIR}"
[ "$SKIP_PHIX" = false ] && echo -e "✓ Will download PhiX174 reference"
[ "$SKIP_HOST" = false ] && echo -e "✓ Will download host genome (${HOST_ORGANISM})"
[ "$SKIP_SILVA" = false ] && echo -e "✓ Will download SILVA rRNA database"
echo -e ""

# Estimate download size and time
echo -e "${YELLOW}Estimated downloads:${NC}"
TOTAL_SIZE="~50 MB"
TOTAL_TIME="~5 minutes"

if [ "$SKIP_PHIX" = false ]; then
    echo -e "  PhiX174: ~5 KB"
fi

if [ "$SKIP_HOST" = false ]; then
    case "$HOST_ORGANISM" in
        human)
            echo -e "  Human genome: ~900 MB compressed, ~3.1 GB uncompressed"
            TOTAL_SIZE="~950 MB"
            TOTAL_TIME="~15-30 minutes"
            ;;
        mouse)
            echo -e "  Mouse genome: ~800 MB compressed, ~2.7 GB uncompressed"
            TOTAL_SIZE="~850 MB"
            TOTAL_TIME="~15-30 minutes"
            ;;
    esac
fi

if [ "$SKIP_SILVA" = false ]; then
    echo -e "  SILVA rRNA: ~100-150 MB compressed, ~350 MB uncompressed"
fi

echo -e "\nTotal: ${TOTAL_SIZE} (${TOTAL_TIME} depending on connection)"
echo -e ""

# Confirm with user
read -p "Continue with setup? (y/N) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo -e "${YELLOW}Setup cancelled${NC}"
    exit 0
fi

echo -e "\n${GREEN}Starting reference database setup...${NC}\n"

# Track success/failure
FAILED=()
SUCCEEDED=()

# Download PhiX174
if [ "$SKIP_PHIX" = false ]; then
    echo -e "${GREEN}[1/3] PhiX174 Reference${NC}"
    echo -e "${GREEN}================================${NC}"
    if bash "${RESOURCES_DIR}/download_phix.sh"; then
        SUCCEEDED+=("PhiX174")
    else
        FAILED+=("PhiX174")
        echo -e "${RED}PhiX174 download failed${NC}"
    fi
    echo -e ""
fi

# Download Host Genome
if [ "$SKIP_HOST" = false ]; then
    echo -e "${GREEN}[2/3] Host Genome Reference${NC}"
    echo -e "${GREEN}================================${NC}"
    if bash "${RESOURCES_DIR}/download_host.sh" "$HOST_ORGANISM"; then
        SUCCEEDED+=("Host genome ($HOST_ORGANISM)")
    else
        FAILED+=("Host genome ($HOST_ORGANISM)")
        echo -e "${RED}Host genome download failed${NC}"
    fi
    echo -e ""
fi

# Download SILVA rRNA
if [ "$SKIP_SILVA" = false ]; then
    echo -e "${GREEN}[3/3] SILVA rRNA Database${NC}"
    echo -e "${GREEN}================================${NC}"
    if bash "${RESOURCES_DIR}/download_silva.sh"; then
        SUCCEEDED+=("SILVA rRNA")
    else
        FAILED+=("SILVA rRNA")
        echo -e "${RED}SILVA download failed${NC}"
    fi
    echo -e ""
fi

# Create version file
echo -e "${GREEN}Creating version tracking file...${NC}"
VERSION_FILE="${RESOURCES_DIR}/VERSIONS.txt"
{
    echo "Lab Virome QC - Reference Database Versions"
    echo "============================================"
    echo ""
    echo "Setup date: $(date)"
    echo ""
    [ -f "${RESOURCES_DIR}/phix174.fasta" ] && echo "PhiX174: NC_001422.1 (NCBI RefSeq)"
    if [ -f "${RESOURCES_DIR}/host_genome.fasta" ]; then
        case "$HOST_ORGANISM" in
            human) echo "Host genome: Homo sapiens GRCh38 (Ensembl release 110)" ;;
            mouse) echo "Host genome: Mus musculus GRCm39 (Ensembl release 110)" ;;
            *) echo "Host genome: Custom ($HOST_ORGANISM)" ;;
        esac
    fi
    [ -f "${RESOURCES_DIR}/silva_rrna.fasta" ] && echo "SILVA rRNA: v138.1 (SSU+LSU, NR99)"
} > "$VERSION_FILE"

# Final summary
echo -e "\n${GREEN}========================================${NC}"
echo -e "${GREEN}Setup Summary${NC}"
echo -e "${GREEN}========================================${NC}"

if [ ${#SUCCEEDED[@]} -gt 0 ]; then
    echo -e "${GREEN}✓ Successfully downloaded:${NC}"
    for item in "${SUCCEEDED[@]}"; do
        echo -e "  - $item"
    done
fi

if [ ${#FAILED[@]} -gt 0 ]; then
    echo -e "\n${RED}✗ Failed to download:${NC}"
    for item in "${FAILED[@]}"; do
        echo -e "  - $item"
    done
fi

echo -e "\n${BLUE}Files in resources/:${NC}"
ls -lh "${RESOURCES_DIR}"/*.fasta 2>/dev/null || echo "No FASTA files found"

# Verification
echo -e "\n${GREEN}Verifying setup...${NC}"
MISSING=()
[ ! -f "${RESOURCES_DIR}/phix174.fasta" ] && MISSING+=("phix174.fasta")
[ ! -f "${RESOURCES_DIR}/host_genome.fasta" ] && MISSING+=("host_genome.fasta")
[ ! -f "${RESOURCES_DIR}/silva_rrna.fasta" ] && MISSING+=("silva_rrna.fasta")

if [ ${#MISSING[@]} -eq 0 ]; then
    echo -e "${GREEN}✓ All required reference files are present!${NC}"
    echo -e "\n${GREEN}You can now run the pipeline:${NC}"
    echo -e "  snakemake --use-conda --cores 8"
else
    echo -e "${YELLOW}⚠ Missing files:${NC}"
    for file in "${MISSING[@]}"; do
        echo -e "  - $file"
    done
    echo -e "\n${YELLOW}Run the individual download scripts to get missing files:${NC}"
    [ ! -f "${RESOURCES_DIR}/phix174.fasta" ] && echo -e "  bash resources/download_phix.sh"
    [ ! -f "${RESOURCES_DIR}/host_genome.fasta" ] && echo -e "  bash resources/download_host.sh <organism>"
    [ ! -f "${RESOURCES_DIR}/silva_rrna.fasta" ] && echo -e "  bash resources/download_silva.sh"
fi

echo -e "\n${GREEN}Setup complete!${NC}"
