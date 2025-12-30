#!/bin/bash
#SBATCH --job-name=hisat2_build_index
#SBATCH --time=4:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=guest

# Build HISAT2 index for Necator americanus reference genome

BASEDIR="/work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data"
REFERENCEDIR="${BASEDIR}/reference"

# Reference genome file
REFERENCE_FASTA="MaSuRCA_config_purged_namericanus_withMito.short.masked.fasta"
INDEX_PREFIX="MaSuRCA_config_purged_namericanus_withMito.short.masked"

echo "==================================================================="
echo "Building HISAT2 index"
echo "Input: ${REFERENCEDIR}/${REFERENCE_FASTA}"
echo "Output prefix: ${REFERENCEDIR}/${INDEX_PREFIX}"
echo "Started at: $(date)"
echo "==================================================================="
printf "\n"

# Check if reference exists
if [ ! -f "${REFERENCEDIR}/${REFERENCE_FASTA}" ]; then
    echo "✗ Error: Reference file not found: ${REFERENCEDIR}/${REFERENCE_FASTA}"
    exit 1
fi

# Load HISAT2 module
module purge
module load hisat2/2.2

# Build HISAT2 index
echo "Building index (this may take a while)..."
hisat2-build \
    -p 8 \
    "${REFERENCEDIR}/${REFERENCE_FASTA}" \
    "${REFERENCEDIR}/${INDEX_PREFIX}"

# Verify index files were created
if ls "${REFERENCEDIR}/${INDEX_PREFIX}".*.ht2 1> /dev/null 2>&1; then
    printf "\n"
    echo "==================================================================="
    echo "✓ HISAT2 index built successfully"
    echo "==================================================================="
    echo "Index files:"
    ls -lh "${REFERENCEDIR}/${INDEX_PREFIX}".*.ht2
    echo "==================================================================="
else
    echo "✗ Error: Index files not created"
    exit 1
fi

echo "Completed at: $(date)"
