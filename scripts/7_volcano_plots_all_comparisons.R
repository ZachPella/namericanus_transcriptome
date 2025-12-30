#!/usr/bin/env Rscript

# ==============================================================================
# Volcano Plots for All Comparisons (matching Adult vs L3 style)
# ==============================================================================

library(ggplot2)
library(dplyr)

BASE_DIR <- "/work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data"
setwd(paste0(BASE_DIR, "/differential_expression"))

cat("================================================================\n")
cat("Creating Volcano Plots for All Comparisons\n")
cat("================================================================\n\n")

# Function to create volcano plot
create_volcano <- function(results_file, comparison_name, title_text) {
    
    cat(paste("Creating volcano plot for", comparison_name, "...\n"))
    
    # Load results
    results <- read.csv(results_file)
    
    # Get comparison parts for labeling
    comparison_parts <- strsplit(comparison_name, "_vs_")[[1]]
    group1 <- comparison_parts[1]
    group2 <- comparison_parts[2]
    
    # Add significance labels
    results <- results %>%
        mutate(
            significant = case_when(
                padj < 0.05 & log2FoldChange > 1 ~ paste0("Upregulated in ", group1),
                padj < 0.05 & log2FoldChange < -1 ~ paste0("Upregulated in ", group2),
                padj < 0.05 ~ "Significant (|FC| < 2)",
                TRUE ~ "Not Significant"
            ),
            neg_log10_padj = -log10(padj)
        )
    
    # Count categories and add to labels
    sig_counts <- results %>%
        filter(!is.na(padj)) %>%
        count(significant)
    
    cat("\nGene categories:\n")
    print(sig_counts)
    cat("\n")
    
    # Create labels with counts
    label_mapping <- setNames(
        paste0(sig_counts$significant, " (n=", sig_counts$n, ")"),
        sig_counts$significant
    )
    
    results$significant_labeled <- label_mapping[results$significant]
    
    # Set up colors with proper syntax
    color_mapping <- list()
    
    for (i in 1:nrow(sig_counts)) {
        cat_name <- sig_counts$significant[i]
        cat_count <- sig_counts$n[i]
        label <- paste0(cat_name, " (n=", cat_count, ")")
        
        if (cat_name == "Not Significant") {
            color_mapping[[label]] <- "gray70"
        } else if (cat_name == "Significant (|FC| < 2)") {
            color_mapping[[label]] <- "#FFA726"
        } else if (cat_name == paste0("Upregulated in ", group1)) {
            color_mapping[[label]] <- "#D32F2F"
        } else if (cat_name == paste0("Upregulated in ", group2)) {
            color_mapping[[label]] <- "#1976D2"
        }
    }
    
    color_mapping <- unlist(color_mapping)
    
    # Create volcano plot
    p <- ggplot(results, aes(x = log2FoldChange, y = neg_log10_padj, color = significant_labeled)) +
        geom_point(alpha = 0.6, size = 1.5) +
        scale_color_manual(values = color_mapping, name = "") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
        labs(
            title = title_text,
            x = paste0("log2 Fold Change (", gsub("_", " ", comparison_name), ")"),
            y = "-log10(adjusted p-value)"
        ) +
        theme_bw(base_size = 14) +
        theme(
            legend.position = "bottom",
            legend.text = element_text(size = 11),
            plot.title = element_text(face = "bold", hjust = 0.5)
        ) +
        guides(color = guide_legend(override.aes = list(size = 3)))
    
    # Save plots
    ggsave(paste0("volcano_plot_", comparison_name, ".pdf"), p, width = 10, height = 8, dpi = 300)
    ggsave(paste0("volcano_plot_", comparison_name, ".png"), p, width = 10, height = 8, dpi = 300)
    
    cat(paste("  ✓ Saved: volcano_plot_", comparison_name, "\n\n"))
}

# Create volcano plots for all three comparisons
create_volcano("DESeq2_results_Female_vs_L3.csv", 
               "Female_vs_L3",
               "Differential Gene Expression: Female vs L3 Hookworms")

create_volcano("DESeq2_results_Male_vs_L3.csv",
               "Male_vs_L3", 
               "Differential Gene Expression: Male vs L3 Hookworms")

create_volcano("DESeq2_results_Female_vs_Male.csv",
               "Female_vs_Male",
               "Differential Gene Expression: Female vs Male Hookworms")

cat("================================================================\n")
cat("✓ All volcano plots created!\n")
cat("================================================================\n")
