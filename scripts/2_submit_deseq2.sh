#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=DESeq2_analysis
#SBATCH --output=DESeq2.%j.out
#SBATCH --error=DESeq2.%j.err
#SBATCH --partition=guest

echo "=================================================================="
echo "DESeq2 Differential Expression Analysis"
echo "=================================================================="
echo "Start time: $(date)"
echo ""

module load R

# Define paths
DIFF_EXP_DIR="/work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data/differential_expression"

cd $DIFF_EXP_DIR

echo "Working directory: $(pwd)"
echo ""

# Verify input file exists
if [ ! -f "gene_counts_all_samples.txt" ]; then
    echo "ERROR: gene_counts_all_samples.txt not found!"
    exit 1
fi

echo "✓ Found input file: gene_counts_all_samples.txt"
GENE_COUNT=$(tail -n +3 gene_counts_all_samples.txt | wc -l)
echo "  Features to analyze: $GENE_COUNT"
echo ""

# Run DESeq2
echo "Running DESeq2 analysis..."
Rscript DESeq2_hookworm_analysis.R

EXIT_CODE=$?

echo ""
echo "=================================================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ DESeq2 analysis completed successfully"
    echo ""
    echo "Output files:"
    ls -lh DESeq2_*.csv *.pdf 2>/dev/null | awk '{print "  "$9, "("$5")"}'
else
    echo "✗ DESeq2 analysis FAILED"
    exit 1
fi
echo "=================================================================="
echo "End time: $(date)"
