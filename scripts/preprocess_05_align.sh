#!/bin/bash
#SBATCH --job-name=hisat2_alignment_all
#SBATCH --time=6-00:00:00
#SBATCH --output=%x_%j_%a.out
#SBATCH --error=%x_%j_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=120G
#SBATCH --array=1-12
#SBATCH --partition=guest

# Align ALL 12 samples to Necator americanus reference genome using HISAT2

BASEDIR="/work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data"
SAMPLE_LIST="${BASEDIR}/sample_list.txt"
READSDIR="${BASEDIR}/trimmed_reads"
REFERENCEDIR="${BASEDIR}/reference"
WORKDIR="${BASEDIR}/sam_files"

# Reference genome (without extension for HISAT2 index)
REFERENCE="MaSuRCA_config_purged_namericanus_withMito.short.masked"

# Create output directory
mkdir -p ${WORKDIR}

# Get sample name
if [ ! -f "${SAMPLE_LIST}" ]; then
    echo "Error: Sample list file not found: ${SAMPLE_LIST}"
    exit 1
fi

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")

if [[ -z "$SAMPLE" ]]; then
    echo "Error: Empty sample name for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

READS1_TRIMMED="${SAMPLE}_R1_trimmed.fastq.gz"
READS2_TRIMMED="${SAMPLE}_R2_trimmed.fastq.gz"

echo "==================================================================="
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Sample: ${SAMPLE}"
echo "Reference: ${REFERENCEDIR}/${REFERENCE}"
echo "Started at: $(date)"
echo "==================================================================="
printf "\n"

# Check if input files exist
if [ ! -f "${READSDIR}/${READS1_TRIMMED}" ] || [ ! -f "${READSDIR}/${READS2_TRIMMED}" ]; then
    echo "✗ Error: Trimmed read files not found for ${SAMPLE}"
    echo "Expected:"
    echo "  ${READSDIR}/${READS1_TRIMMED}"
    echo "  ${READSDIR}/${READS2_TRIMMED}"
    exit 1
fi

# Check if HISAT2 index exists (check for .1.ht2 file)
if [ ! -f "${REFERENCEDIR}/${REFERENCE}.1.ht2" ]; then
    echo "✗ Error: HISAT2 index not found. Please build index first."
    echo "Expected: ${REFERENCEDIR}/${REFERENCE}.*.ht2"
    exit 1
fi

# Load HISAT2 module
module purge
module load hisat2/2.2

# Run HISAT2 alignment
echo "Running HISAT2 alignment..."
hisat2 \
    -p 16 \
    --dta \
    -x "${REFERENCEDIR}/${REFERENCE}" \
    -1 "${READSDIR}/${READS1_TRIMMED}" \
    -2 "${READSDIR}/${READS2_TRIMMED}" \
    -S "${WORKDIR}/${SAMPLE}.sam" \
    --summary-file "${WORKDIR}/${SAMPLE}.hisat2.summary.txt" \
    2> "${WORKDIR}/${SAMPLE}.hisat2.log"

# Verify output file was created
if [[ -f "${WORKDIR}/${SAMPLE}.sam" && -s "${WORKDIR}/${SAMPLE}.sam" ]]; then
    printf "\n"
    echo "==================================================================="
    echo "✓ HISAT2 alignment completed successfully for ${SAMPLE}"
    echo "==================================================================="
    echo "Output SAM file:"
    ls -lh "${WORKDIR}/${SAMPLE}.sam"
    echo ""
    echo "Alignment summary:"
    cat "${WORKDIR}/${SAMPLE}.hisat2.summary.txt"
    echo "==================================================================="
else
    echo "✗ Error: SAM file not created or is empty for ${SAMPLE}"
    exit 1
fi

echo "Completed at: $(date)"
