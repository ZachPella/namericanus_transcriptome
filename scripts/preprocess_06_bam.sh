#!/bin/bash
#SBATCH --job-name=sam_to_bam_all
#SBATCH --time=4-00:00:00
#SBATCH --output=%x_%j_%a.out
#SBATCH --error=%x_%j_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=80G
#SBATCH --array=1-12
#SBATCH --partition=guest

# Convert SAM to BAM, sort, filter, and generate statistics for ALL 12 samples
# Works with SAM files from HISAT2 (or any other aligner)

BASEDIR="/work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data"
SAMDIR="${BASEDIR}/sam_files"
WORKDIR="${BASEDIR}/bam_files"
REFERENCEDIR="${BASEDIR}/reference"
SAMPLE_LIST="${BASEDIR}/sample_list.txt"

# Reference file
REFERENCE="MaSuRCA_config_purged_namericanus_withMito.short.masked.fasta"
TARGETS="${REFERENCEDIR}/${REFERENCE}.bed"

# Create working directory
mkdir -p "${WORKDIR}"

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

SAM_INPUT="${SAMDIR}/${SAMPLE}.sam"
BAM_UNSORTED="${WORKDIR}/${SAMPLE}.bam"
BAM_SORTED="${WORKDIR}/${SAMPLE}.sorted.bam"
BAM_MAPPED="${WORKDIR}/${SAMPLE}.sorted.mapped.bam"

echo "==================================================================="
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Sample: ${SAMPLE}"
echo "Input SAM: ${SAM_INPUT}"
echo "Started at: $(date)"
echo "==================================================================="
printf "\n"

# Check if SAM file exists
if [ ! -f "${SAM_INPUT}" ]; then
    echo "✗ Error: SAM file not found: ${SAM_INPUT}"
    exit 1
fi

# Load samtools
module purge
module load samtools/1.20

# Convert SAM to BAM
echo "Step 1: Converting SAM to unsorted BAM..."
samtools view -@ 2 -b "${SAM_INPUT}" > "${BAM_UNSORTED}"

# Sort BAM file
echo "Step 2: Sorting BAM file..."
samtools sort -@ 2 -m 16G "${BAM_UNSORTED}" -o "${BAM_SORTED}" -T "${WORKDIR}/${SAMPLE}.reads.tmp"

# Remove intermediate files if successful
if [ -s "${BAM_SORTED}" ]; then
    rm "${SAM_INPUT}"
    rm "${BAM_UNSORTED}"
    echo "  ✓ SAM-to-BAM conversion successful. Intermediate files cleaned."
else
    echo "  ✗ Error: Sorted BAM file is empty or conversion failed"
    exit 1
fi

# Index sorted BAM
echo "Step 3: Indexing sorted BAM..."
samtools index "${BAM_SORTED}"

# Generate statistics
echo "Step 4: Generating alignment statistics..."

# Flagstat
echo "  - Running flagstat..."
samtools flagstat "${BAM_SORTED}" > "${WORKDIR}/flagstat.${SAMPLE}.out"

# Stats with coverage thresholds
if [ -f "${TARGETS}" ]; then
    echo "  - Running stats with coverage thresholds..."
    for cov in 5 10; do
        samtools stats -t "${TARGETS}" --cov-threshold ${cov} "${BAM_SORTED}" > "${WORKDIR}/stats.${cov}x.${SAMPLE}.out"
    done
else
    echo "  - Running general stats (no BED file)..."
    samtools stats "${BAM_SORTED}" > "${WORKDIR}/stats.general.${SAMPLE}.out"
fi

# Depth of coverage
echo "  - Calculating depth of coverage..."
samtools depth -a "${BAM_SORTED}" > "${WORKDIR}/${SAMPLE}.depth"
AVGDOC=$(awk '{ total += $3; count++ } END { print total/count }' "${WORKDIR}/${SAMPLE}.depth")
echo "Average depth of coverage: ${AVGDOC}" > "${WORKDIR}/averageDOC.${SAMPLE}.out"

# Coverage by contig
echo "  - Calculating coverage by contig..."
samtools coverage -o "${WORKDIR}/coverage.${SAMPLE}.out" "${BAM_SORTED}"
samtools coverage --plot-depth -o "${WORKDIR}/hist.coverage.${SAMPLE}.out" "${BAM_SORTED}"

# Create mapped-only BAM
echo "Step 5: Creating mapped-reads-only BAM..."
samtools view -@ 2 -b -F 4 "${BAM_SORTED}" > "${BAM_MAPPED}"
samtools index "${BAM_MAPPED}"

# Final verification
if [[ -f "${BAM_MAPPED}" && -f "${BAM_MAPPED}.bai" ]]; then
    printf "\n"
    echo "==================================================================="
    echo "✓ Processing completed successfully for ${SAMPLE}"
    echo "==================================================================="
    echo "Output files:"
    echo "  Main BAM: ${BAM_SORTED}"
    echo "  Mapped BAM: ${BAM_MAPPED}"
    echo ""
    echo "Statistics files:"
    ls -lh "${WORKDIR}"/*${SAMPLE}*
    echo "==================================================================="
else
    echo "✗ Error: Final BAM files not created for ${SAMPLE}"
    exit 1
fi

echo "Completed at: $(date)"
