#!/usr/bin/env Rscript

# ==============================================================================
# Volcano Plot: Adult vs L3
# ==============================================================================

library(ggplot2)
library(dplyr)

BASE_DIR <- "/work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data"
OUTPUT_DIR <- paste0(BASE_DIR, "/differential_expression_total_adults_vs_L3")

setwd(OUTPUT_DIR)

cat("================================================================\n")
cat("Creating Volcano Plot: Adult vs L3\n")
cat("================================================================\n\n")

# Load DESeq2 results
results <- read.csv("DESeq2_results_Adult_vs_L3.csv")

cat(paste("Total genes:", nrow(results), "\n"))

# Add significance labels
results <- results %>%
    mutate(
        significant = case_when(
            padj < 0.05 & log2FoldChange > 1 ~ "Upregulated in Adult",
            padj < 0.05 & log2FoldChange < -1 ~ "Upregulated in L3",
            padj < 0.05 ~ "Significant (|FC| < 2)",
            TRUE ~ "Not Significant"
        ),
        neg_log10_padj = -log10(padj)
    )

# Count each category
sig_counts <- results %>%
    filter(!is.na(padj)) %>%
    count(significant)

cat("\nGene categories:\n")
print(sig_counts)
cat("\n")

# Create volcano plot
p <- ggplot(results, aes(x = log2FoldChange, y = neg_log10_padj, color = significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(
        values = c(
            "Upregulated in Adult" = "#D32F2F",
            "Upregulated in L3" = "#1976D2",
            "Significant (|FC| < 2)" = "#FFA726",
            "Not Significant" = "gray70"
        ),
        name = ""
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
    labs(
        title = "Differential Gene Expression: Adult vs L3 Hookworms",
        subtitle = paste0("N. americanus transcriptome (n=", nrow(results), " genes)"),
        x = "log2 Fold Change (Adult / L3)",
        y = "-log10(adjusted p-value)"
    ) +
    theme_bw(base_size = 14) +
    theme(
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))

# Save high-resolution plot
ggsave("volcano_plot_Adult_vs_L3.pdf", p, width = 10, height = 8, dpi = 300)
ggsave("volcano_plot_Adult_vs_L3.png", p, width = 10, height = 8, dpi = 300)

cat("================================================================\n")
cat("✓ Volcano plot saved!\n")
cat("  - volcano_plot_Adult_vs_L3.pdf\n")
cat("  - volcano_plot_Adult_vs_L3.png\n")
cat("================================================================\n")
