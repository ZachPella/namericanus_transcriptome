#!/bin/bash
#SBATCH --job-name=concatenate_lanes_all
#SBATCH --time=3-00:00:00
#SBATCH --output=%x_%j_%a.out
#SBATCH --error=%x_%j_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --array=1-12
#SBATCH --partition=guest

# Concatenate lanes for ALL samples (old + new)
# OLD data: Already single lane, just copy/rename
# NEW data: Concatenate L1 + L2 lanes

BASE_DIR="/work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data"
RAW_FASTQ_DIR="${BASE_DIR}/raw_fastq"
OUTPUT_DIR="${BASE_DIR}/concatenated_reads"
SAMPLE_LIST="${BASE_DIR}/sample_list.txt"

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Check if sample list exists
if [ ! -f "$SAMPLE_LIST" ]; then
    echo "Error: Sample list file not found: $SAMPLE_LIST"
    exit 1
fi

# Get sample name from list
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

if [[ -z "$SAMPLE_NAME" ]]; then
    echo "Error: Empty sample name for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Processing sample: ${SAMPLE_NAME}"
echo "Started at: $(date)"
printf "\n"

# Determine if this is OLD or NEW data based on sample name
if [[ "$SAMPLE_NAME" == Na-* ]]; then
    # OLD DATA: Single lane, just copy/rename
    echo "Processing OLD data (single lane)..."
    
    R1_FILE="${RAW_FASTQ_DIR}/${SAMPLE_NAME}_R1_001.fastq.gz"
    R2_FILE="${RAW_FASTQ_DIR}/${SAMPLE_NAME}_R2_001.fastq.gz"
    
    if [[ -f "$R1_FILE" && -f "$R2_FILE" ]]; then
        echo "Found files:"
        echo "  R1: $R1_FILE"
        echo "  R2: $R2_FILE"
        
        # Copy (or symlink) to concatenated directory with consistent naming
        echo "Copying to concatenated_reads directory..."
        cp "$R1_FILE" "${OUTPUT_DIR}/${SAMPLE_NAME}_R1_merged.fastq.gz"
        cp "$R2_FILE" "${OUTPUT_DIR}/${SAMPLE_NAME}_R2_merged.fastq.gz"
        
        echo "✓ Successfully processed ${SAMPLE_NAME} (old data)"
    else
        echo "✗ Missing files for ${SAMPLE_NAME}"
        echo "Looking for: $R1_FILE"
        exit 1
    fi
    
else
    # NEW DATA: Dual lanes (L1 + L2), need concatenation
    echo "Processing NEW data (dual lanes L1 + L2)..."
    
    L1_R1="${RAW_FASTQ_DIR}/${SAMPLE_NAME}_L1_R1.fastq.gz"
    L1_R2="${RAW_FASTQ_DIR}/${SAMPLE_NAME}_L1_R2.fastq.gz"
    L2_R1="${RAW_FASTQ_DIR}/${SAMPLE_NAME}_L2_R1.fastq.gz"
    L2_R2="${RAW_FASTQ_DIR}/${SAMPLE_NAME}_L2_R2.fastq.gz"
    
    if [[ -f "$L1_R1" && -f "$L1_R2" && -f "$L2_R1" && -f "$L2_R2" ]]; then
        echo "Found all required files:"
        echo "  L1_R1: $L1_R1"
        echo "  L1_R2: $L1_R2"
        echo "  L2_R1: $L2_R1"
        echo "  L2_R2: $L2_R2"
        
        echo "Concatenating R1 files..."
        cat "$L1_R1" "$L2_R1" > "${OUTPUT_DIR}/${SAMPLE_NAME}_R1_merged.fastq.gz"
        
        echo "Concatenating R2 files..."
        cat "$L1_R2" "$L2_R2" > "${OUTPUT_DIR}/${SAMPLE_NAME}_R2_merged.fastq.gz"
        
        echo "✓ Successfully concatenated ${SAMPLE_NAME} (new data)"
    else
        echo "✗ Missing required files for ${SAMPLE_NAME}"
        echo "Looking in: $RAW_FASTQ_DIR"
        ls -la "${RAW_FASTQ_DIR}/${SAMPLE_NAME}"* 2>/dev/null || echo "No matching files found"
        exit 1
    fi
fi

# Verify output files were created
if [[ -f "${OUTPUT_DIR}/${SAMPLE_NAME}_R1_merged.fastq.gz" && \
      -f "${OUTPUT_DIR}/${SAMPLE_NAME}_R2_merged.fastq.gz" ]]; then
    echo ""
    echo "=== Output verification ==="
    ls -lh "${OUTPUT_DIR}/${SAMPLE_NAME}"*merged.fastq.gz
    echo "✓ Output files created successfully"
else
    echo "✗ Error: Output files not created for ${SAMPLE_NAME}"
    exit 1
fi

echo ""
echo "Completed at: $(date)"
