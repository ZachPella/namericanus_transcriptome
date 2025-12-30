#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=88G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=interproscan
#SBATCH --output=interproscan.%j.out
#SBATCH --error=interproscan.%j.err
#SBATCH --partition=guest

echo "================================================================"
echo "Running InterProScan on BRAKER3 proteins"
echo "================================================================"
echo "Start time: $(date)"
echo ""

module load interproscan/5.69
module load java/11

# Directories
BASE_DIR="/work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data"
PROTEIN_FILE="${BASE_DIR}/braker_annotation/braker_output/Augustus/augustus.hints.aa"
OUTPUT_DIR="${BASE_DIR}/functional_annotation/interproscan"
OUTPUT_BASE="namericanus_interproscan"

mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

echo "Input: $PROTEIN_FILE"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Check input
if [ ! -f "$PROTEIN_FILE" ]; then
    echo "ERROR: Protein file not found!"
    exit 1
fi

# Count proteins
PROTEIN_COUNT=$(grep -c ">" $PROTEIN_FILE)
echo "Proteins to annotate: $PROTEIN_COUNT"
echo ""

# Run InterProScan
# -appl: applications to run (using most informative ones)
# -f: output formats
# -goterms: include GO term annotations
# -pa: include pathway annotations
# -cpu: number of threads

echo "Running InterProScan (this will take 24-48 hours)..."
interproscan.sh \
    -i $PROTEIN_FILE \
    -f tsv,gff3 \
    -goterms \
    -pa \
    -cpu 16 \
    -appl Pfam,SMART,SUPERFAMILY,Gene3D,PRINTS,ProSiteProfiles,Hamap

EXIT_CODE=$?

echo ""
echo "================================================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ InterProScan completed successfully"
    echo ""
    echo "Output files:"
    ls -lh ${OUTPUT_BASE}*
    echo ""
    
    # Count annotations
    if [ -f "${OUTPUT_BASE}.tsv" ]; then
        ANNOTATED=$(cut -f1 ${OUTPUT_BASE}.tsv | sort -u | wc -l)
        TOTAL_ANNOT=$(wc -l < ${OUTPUT_BASE}.tsv)
        echo "Statistics:"
        echo "  Unique proteins annotated: $ANNOTATED / $PROTEIN_COUNT"
        echo "  Total annotations: $TOTAL_ANNOT"
    fi
else
    echo "✗ InterProScan FAILED"
    exit 1
fi

echo "================================================================"
echo "End time: $(date)"
