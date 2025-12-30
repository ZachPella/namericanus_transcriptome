#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=DESeq2_adults_L3
#SBATCH --output=DESeq2_adults_L3.%j.out
#SBATCH --error=DESeq2_adults_L3.%j.err

echo "=================================================================="
echo "DESeq2: Adult (combined) vs L3 Analysis"
echo "=================================================================="
echo "Start time: $(date)"

module load R

BASE_DIR="/work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data"
INPUT_DIR="${BASE_DIR}/differential_expression"
OUTPUT_DIR="${BASE_DIR}/differential_expression_total_adults_vs_L3"

mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

echo "Input: ${INPUT_DIR}/gene_counts_all_samples.txt"
echo "Output: ${OUTPUT_DIR}"
echo ""

Rscript - <<'EOF'
library(DESeq2)
library(dplyr)

BASE_DIR <- "/work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data"
INPUT_DIR <- paste0(BASE_DIR, "/differential_expression")
OUTPUT_DIR <- paste0(BASE_DIR, "/differential_expression_total_adults_vs_L3")

setwd(OUTPUT_DIR)

# Load featureCounts output from original pipeline
counts_raw <- read.table(paste0(INPUT_DIR, "/gene_counts_all_samples.txt"),
                         header = TRUE,
                         row.names = 1,
                         skip = 1)

# Remove non-count columns
counts <- counts_raw[, -c(1:5)]

# Simplify column names
colnames(counts) <- gsub(".*/", "", colnames(counts))
colnames(counts) <- gsub(".sorted.bam", "", colnames(counts))

cat("Sample names:\n")
print(colnames(counts))
cat("\n")

# Create sample metadata with Adult group
sample_names <- colnames(counts)
coldata <- data.frame(
    sample = sample_names,
    stage = c("Female", "Female", "Female", "Female",
              "Male", "Male", "Male", "Male",
              "L3", "L3", "L3", "L3"),
    stringsAsFactors = FALSE
)

# Combine Female and Male into Adult
coldata$group <- ifelse(coldata$stage == "L3", "L3", "Adult")
rownames(coldata) <- sample_names

cat("Sample metadata:\n")
print(coldata)
cat("\n")

# ============================================================================
# Adult vs L3 Comparison
# ============================================================================

cat("Running Adult vs L3 analysis...\n\n")

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ group)

# Filter low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

cat(paste("Genes after filtering:", nrow(dds), "\n\n"))

# Set L3 as reference BEFORE running DESeq
dds$group <- relevel(dds$group, ref = "L3")

# Run DESeq2
dds <- DESeq(dds)

# Get results: Adult vs L3
res <- results(dds)
res <- as.data.frame(res)
res$gene_id <- rownames(res)

# Reorder columns
res <- res[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

# Filter significant genes
res_sig <- res[res$padj < 0.05 & !is.na(res$padj), ]

# Save results
write.csv(res, "DESeq2_results_Adult_vs_L3.csv", row.names = FALSE)
write.csv(res_sig, "DESeq2_significant_Adult_vs_L3.csv", row.names = FALSE)

# Summary statistics
cat("================================================================\n")
cat("RESULTS SUMMARY\n")
cat("================================================================\n\n")

cat(paste("Total genes tested:", nrow(res), "\n"))
cat(paste("Significant genes (padj < 0.05):", nrow(res_sig), "\n"))
cat(paste("  Upregulated in Adult:", sum(res_sig$log2FoldChange > 0), "\n"))
cat(paste("  Downregulated in Adult:", sum(res_sig$log2FoldChange < 0), "\n\n"))

# Top upregulated genes
cat("Top 10 upregulated in Adult:\n")
top_up <- head(res_sig[order(-res_sig$log2FoldChange), c("gene_id", "log2FoldChange", "padj")], 10)
print(top_up)
cat("\n")

# Top downregulated genes
cat("Top 10 downregulated in Adult (upregulated in L3):\n")
top_down <- head(res_sig[order(res_sig$log2FoldChange), c("gene_id", "log2FoldChange", "padj")], 10)
print(top_down)
cat("\n")

cat("================================================================\n")
cat("✓ Analysis complete!\n")
cat("================================================================\n")

EOF

EXIT_CODE=$?

echo ""
echo "=================================================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ DESeq2 analysis completed successfully"
    echo ""
    echo "Output files in: ${OUTPUT_DIR}"
    ls -lh ${OUTPUT_DIR}/*.csv
else
    echo "✗ DESeq2 analysis FAILED"
    exit 1
fi
echo "=================================================================="
echo "End time: $(date)"
