#!/bin/bash
#SBATCH --job-name=fastp_trimming_all
#SBATCH --time=1-06:00:00
#SBATCH --output=%x_%j_%a.out
#SBATCH --error=%x_%j_%a.err
#SBATCH --ntasks=1
#SBATCH --mem=45G
#SBATCH --ntasks-per-node=8
#SBATCH --array=1-12
#SBATCH --partition=guest

# Quality control and adapter trimming for ALL 12 samples

BASEDIR="/work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data"
READSDIR="${BASEDIR}/concatenated_reads"
WORKDIR="${BASEDIR}/trimmed_reads"
SAMPLE_LIST="${BASEDIR}/sample_list.txt"

# Create output directory
mkdir -p "${WORKDIR}"

# Get sample name
if [ ! -f "${SAMPLE_LIST}" ]; then
    echo "Error: Sample list file not found: ${SAMPLE_LIST}"
    exit 1
fi

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")

if [[ -z "$SAMPLE" ]]; then
    echo "Error: Sample name is empty for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

READS1="${SAMPLE}_R1_merged.fastq.gz"
READS2="${SAMPLE}_R2_merged.fastq.gz"
READS1_TRIMMED="${SAMPLE}_R1_trimmed.fastq.gz"
READS2_TRIMMED="${SAMPLE}_R2_trimmed.fastq.gz"

echo "==================================================================="
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Processing Sample: ${SAMPLE}"
echo "Input R1: ${READSDIR}/${READS1}"
echo "Input R2: ${READSDIR}/${READS2}"
echo "Started at: $(date)"
echo "==================================================================="
printf "\n"

# Verify input files exist
if [ ! -f "${READSDIR}/${READS1}" ] || [ ! -f "${READSDIR}/${READS2}" ]; then
    echo "✗ Error: Concatenated input files not found for ${SAMPLE}"
    echo "Expected:"
    echo "  ${READSDIR}/${READS1}"
    echo "  ${READSDIR}/${READS2}"
    exit 1
fi

# Load fastp module
module purge
module load fastp/0.23

# Change to working directory
cd "${WORKDIR}"

# Run fastp
echo "Running fastp..."
fastp \
    --in1 "${READSDIR}/${READS1}" \
    --in2 "${READSDIR}/${READS2}" \
    --out1 "${READS1_TRIMMED}" \
    --out2 "${READS2_TRIMMED}" \
    -l 50 \
    -h "${SAMPLE}.fastp.html" \
    -j "${SAMPLE}.fastp.json" \
    --thread 8

# Verify output files were created
if [[ -f "${WORKDIR}/${READS1_TRIMMED}" && -f "${WORKDIR}/${READS2_TRIMMED}" ]]; then
    printf "\n"
    echo "==================================================================="
    echo "✓ fastp completed successfully for ${SAMPLE}"
    echo "==================================================================="
    echo "Output files:"
    ls -lh "${WORKDIR}/${SAMPLE}"*trimmed.fastq.gz
    echo ""
    echo "QC reports:"
    ls -lh "${WORKDIR}/${SAMPLE}".fastp.*
    echo "==================================================================="
else
    echo "✗ Error: fastp output files not created for ${SAMPLE}"
    exit 1
fi

echo "Completed at: $(date)"
