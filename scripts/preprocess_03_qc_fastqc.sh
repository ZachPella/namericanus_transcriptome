#!/bin/bash
#SBATCH --job-name=fastqc_reports_all
#SBATCH --time=1-06:00:00
#SBATCH --output=%x_%j_%a.out
#SBATCH --error=%x_%j_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
#SBATCH --array=1-12
#SBATCH --partition=guest

# Generate FastQC quality reports for ALL 12 samples

BASEDIR="/work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data"
READSDIR="${BASEDIR}/trimmed_reads"
QCDIR="${BASEDIR}/fastqc_reports"
SAMPLE_LIST="${BASEDIR}/sample_list.txt"

# Create QC directory
mkdir -p "${QCDIR}"

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
echo "Started at: $(date)"
echo "==================================================================="
printf "\n"

# Check if input files exist
if [[ ! -f "${READSDIR}/${READS1_TRIMMED}" || ! -f "${READSDIR}/${READS2_TRIMMED}" ]]; then
    echo "✗ ERROR: Trimmed files not found for ${SAMPLE}"
    echo "Expected:"
    echo "  ${READSDIR}/${READS1_TRIMMED}"
    echo "  ${READSDIR}/${READS2_TRIMMED}"
    exit 1
fi

# Load FastQC module
module purge
module load fastqc/0.12

# Run FastQC on both R1 and R2
echo "Running FastQC on R1..."
fastqc --threads 4 --outdir="${QCDIR}" "${READSDIR}/${READS1_TRIMMED}"

echo "Running FastQC on R2..."
fastqc --threads 4 --outdir="${QCDIR}" "${READSDIR}/${READS2_TRIMMED}"

# Verify output files were created
R1_HTML="${QCDIR}/${SAMPLE}_R1_trimmed_fastqc.html"
R2_HTML="${QCDIR}/${SAMPLE}_R2_trimmed_fastqc.html"

if [[ -f "${R1_HTML}" && -f "${R2_HTML}" ]]; then
    printf "\n"
    echo "==================================================================="
    echo "✓ FastQC completed successfully for ${SAMPLE}"
    echo "==================================================================="
    echo "Output files:"
    ls -lh "${QCDIR}/${SAMPLE}"*_fastqc.*
    echo "==================================================================="
else
    echo "✗ Error: FastQC output files not created for ${SAMPLE}"
    exit 1
fi

echo "Completed at: $(date)"
