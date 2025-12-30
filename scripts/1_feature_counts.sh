#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=featureCounts_all_samples
#SBATCH --output=featureCounts.%j.out
#SBATCH --error=featureCounts.%j.err
#SBATCH --partition=guest

echo "=================================================================="
echo "Running featureCounts on all 12 hookworm samples"
echo "=================================================================="
echo "Start time: $(date)"

module load subread

# Directories
BAM_DIR="/work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data/bam_files"
GTF_FILE="/work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data/braker_annotation/braker_output/Augustus/augustus.hints.gtf"
OUT_DIR="/work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data/differential_expression"

mkdir -p $OUT_DIR
cd $OUT_DIR

echo "GTF file: $GTF_FILE"
echo "Output directory: $OUT_DIR"
echo ""

# Verify GTF exists
if [ ! -f "$GTF_FILE" ]; then
    echo "ERROR: GTF file not found at $GTF_FILE"
    echo "Checking what's available in braker_output:"
    ls -lh /work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data/braker_annotation/braker_output/*.gtf
    ls -lh /work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data/braker_annotation/braker_output/Augustus/*.gtf
    exit 1
fi

echo "✓ GTF file found: $(ls -lh $GTF_FILE | awk '{print $5}')"
echo ""

# List all BAM files
echo "BAM files to process:"
ls -1 $BAM_DIR/*.sorted.bam | wc -l
echo " BAM files found"
echo ""

# Verify BAM files exist
if [ ! -d "$BAM_DIR" ]; then
    echo "ERROR: BAM directory not found!"
    exit 1
fi

BAM_COUNT=$(ls -1 $BAM_DIR/*.sorted.bam 2>/dev/null | wc -l)
if [ $BAM_COUNT -eq 0 ]; then
    echo "ERROR: No BAM files found in $BAM_DIR"
    exit 1
fi

echo "Found $BAM_COUNT BAM files"
echo ""

# Run featureCounts with proper settings for RNA-seq
echo "Running featureCounts..."
echo ""

featureCounts \
    -p \
    -B \
    -C \
    --primary \
    -T 16 \
    -t exon \
    -g gene_id \
    -a $GTF_FILE \
    -o gene_counts_all_samples.txt \
    $BAM_DIR/Na-2-female-mRNA_S2_L005.sorted.bam \
    $BAM_DIR/2femaleB.sorted.bam \
    $BAM_DIR/2femaleC.sorted.bam \
    $BAM_DIR/2femaleD.sorted.bam \
    $BAM_DIR/Na-3-male-mRNA_S1_L005.sorted.bam \
    $BAM_DIR/3maleB.sorted.bam \
    $BAM_DIR/3maleC.sorted.bam \
    $BAM_DIR/3maleD.sorted.bam \
    $BAM_DIR/Na-L3-A_AiD-mRNA_S3_L005.sorted.bam \
    $BAM_DIR/L3E.sorted.bam \
    $BAM_DIR/L3F.sorted.bam \
    $BAM_DIR/L3G.sorted.bam

EXIT_CODE=$?

echo ""
echo "=================================================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ featureCounts completed successfully"
    echo ""
    echo "Output file: $(pwd)/gene_counts_all_samples.txt"
    echo "Summary file: $(pwd)/gene_counts_all_samples.txt.summary"
    echo ""
    
    # Show quick stats
    if [ -f "gene_counts_all_samples.txt.summary" ]; then
        echo "Alignment summary:"
        cat gene_counts_all_samples.txt.summary
        echo ""
    fi
    
    # Count features
    if [ -f "gene_counts_all_samples.txt" ]; then
        FEATURE_COUNT=$(tail -n +3 gene_counts_all_samples.txt | wc -l)
        echo "Total features quantified: $FEATURE_COUNT"
    fi
    
    echo ""
    echo "Ready for DESeq2 analysis!"
    echo "Next step: sbatch run_DESeq2.sh"
else
    echo "✗ featureCounts FAILED with exit code $EXIT_CODE"
    exit 1
fi

echo "=================================================================="
echo "End time: $(date)"
