#!/bin/bash
#SBATCH --time=7-00:00:00
#SBATCH --mem=150G
#SBATCH --job-name=braker3_singularity
#SBATCH --output=braker3_singularity.%j.out
#SBATCH --error=braker3_singularity.%j.err
#SBATCH --cpus-per-task=16
#SBATCH --partition=guest

echo "==================================================================="
echo "=== BRAKER3 Annotation with Singularity Container ==="
echo "==================================================================="
echo "Start time: $(date)"
echo "Host: $(hostname)"
printf "\n"

# Directory setup
BASE_DIR="/work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data"
BRAKER_DIR="${BASE_DIR}/braker_annotation"
BAM_DIR="${BASE_DIR}/bam_files"
CONTAINER="${BASE_DIR}/containers/braker3.sif"
AUGUSTUS_CONFIG="${BRAKER_DIR}/augustus_config"

# Input files
GENOME="${BASE_DIR}/reference/MaSuRCA_config_purged_namericanus_withMito.short.masked.fasta"
PROTEIN_DB="${BASE_DIR}/protein_for_braker3/necator_americanus.PRJNA72135.WBPS19.protein.cleaned.fa"

echo "=== Configuration ==="
echo "Container: $(basename $CONTAINER)"
echo "Genome: $(basename $GENOME)"
echo "Proteins: $(basename $PROTEIN_DB)"
echo "BAM directory: $BAM_DIR"
echo "Augustus config: $AUGUSTUS_CONFIG"
printf "\n"

# Verify basic inputs
echo "Verifying inputs..."
for file in "$CONTAINER" "$GENOME" "$PROTEIN_DB" "${HOME}/.gm_key"; do
    if [ ! -f "$file" ]; then
        echo "✗ Error: File not found: $file"
        exit 1
    fi
done

# Build comma-separated list of BAM files
echo "Finding BAM files..."
BAM_FILES=$(find "$BAM_DIR" -name "*.sorted.bam" | tr '\n' ',' | sed 's/,$//')

if [ -z "$BAM_FILES" ]; then
    echo "✗ Error: No sorted BAM files found in $BAM_DIR"
    exit 1
fi

# Count BAM files
BAM_COUNT=$(echo "$BAM_FILES" | tr ',' '\n' | wc -l)
echo "  ✓ Found $BAM_COUNT BAM files"

# Display first few BAM files
echo "  BAM files to use:"
echo "$BAM_FILES" | tr ',' '\n' | head -5 | sed 's/^/    /'
if [ $BAM_COUNT -gt 5 ]; then
    echo "    ... and $((BAM_COUNT - 5)) more"
fi
printf "\n"

# Verify Augustus config
if [ ! -f "$AUGUSTUS_CONFIG/species/generic/generic_parameters.cfg" ]; then
    echo "✗ Error: Augustus config incomplete - missing generic template"
    exit 1
fi

echo "  ✓ All inputs verified"
echo "  ✓ Augustus config complete"
printf "\n"

# Load singularity
module load singularity

# Clean up any previous output directory
if [ -d "${BRAKER_DIR}/braker_output" ]; then
    echo "Removing previous BRAKER output directory..."
    rm -rf "${BRAKER_DIR}/braker_output"
fi

echo "==================================================================="
echo "Running BRAKER3 in Singularity container..."
echo "Expected runtime: 4-8 hours"
echo "==================================================================="
printf "\n"

# Run BRAKER3
singularity exec \
    --env AUGUSTUS_CONFIG_PATH="${AUGUSTUS_CONFIG}" \
    -B "${BASE_DIR}:${BASE_DIR}" \
    "$CONTAINER" \
    braker.pl \
        --AUGUSTUS_CONFIG_PATH="${AUGUSTUS_CONFIG}" \
        --genome="$GENOME" \
        --prot_seq="$PROTEIN_DB" \
        --bam="$BAM_FILES" \
        --species=namericanus_hisat2_all_12_samples_Dec2025 \
        --workingdir="${BRAKER_DIR}/braker_output" \
        --softmasking \
        --gff3 \
        --verbosity=3 \
        --threads=16

BRAKER_EXIT_CODE=$?

printf "\n"
echo "==================================================================="
if [ $BRAKER_EXIT_CODE -eq 0 ]; then
    echo "✓✓✓ BRAKER3 ANNOTATION COMPLETED SUCCESSFULLY! ✓✓✓"
    echo "==================================================================="
    printf "\n"

    echo "=== Output Files ==="
    echo "Working directory: ${BRAKER_DIR}/braker_output"
    printf "\n"

    # List key output files
    if [ -d "${BRAKER_DIR}/braker_output" ]; then
        echo "Key files:"
        ls -lh "${BRAKER_DIR}/braker_output"/*.gtf 2>/dev/null
        ls -lh "${BRAKER_DIR}/braker_output"/*.gff3 2>/dev/null
        ls -lh "${BRAKER_DIR}/braker_output"/*.aa 2>/dev/null
        printf "\n"
    fi

    # Count genes if GTF exists
    if [ -f "${BRAKER_DIR}/braker_output/augustus.hints.gtf" ]; then
        GENE_COUNT=$(grep -c $'\tgene\t' "${BRAKER_DIR}/braker_output/augustus.hints.gtf" || echo "0")
        TRANSCRIPT_COUNT=$(grep -c $'\ttranscript\t' "${BRAKER_DIR}/braker_output/augustus.hints.gtf" || echo "0")
        echo "=== Summary Statistics ==="
        echo "  Predicted genes:       $GENE_COUNT"
        echo "  Predicted transcripts: $TRANSCRIPT_COUNT"
        printf "\n"
    fi

    echo "=== Next Steps ==="
    echo "1. Run featureCounts on all 12 BAM files using:"
    echo "   ${BRAKER_DIR}/braker_output/augustus.hints.gtf"
    echo "2. Proceed with DESeq2 differential expression analysis"
    printf "\n"

else
    echo "✗✗✗ BRAKER3 FAILED ✗✗✗"
    echo "==================================================================="
    echo "Exit code: $BRAKER_EXIT_CODE"
    printf "\n"
    echo "Check error details in:"
    echo "  - braker3_singularity.*.err"
    echo "  - ${BRAKER_DIR}/braker_output/braker.log"
    echo "  - ${BRAKER_DIR}/braker_output/errors/"
    exit 1
fi

echo "==================================================================="
echo "Pipeline completed at: $(date)"
echo "==================================================================="
