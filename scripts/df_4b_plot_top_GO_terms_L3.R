#!/usr/bin/env Rscript

# ==============================================================================
# Top 10 GO Terms Bar Charts: L3-UPREGULATED
# ==============================================================================

library(ggplot2)
library(dplyr)

BASE_DIR <- "/work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data"
OUTPUT_DIR <- paste0(BASE_DIR, "/differential_expression_total_adults_vs_L3")

setwd(OUTPUT_DIR)

cat("================================================================\n")
cat("Creating Top 10 GO Term Bar Charts - L3 Upregulated\n")
cat("================================================================\n\n")

# Function to create bar chart for top 10 GO terms
plot_top_go_terms <- function(go_file, ontology_name, output_name) {

    cat(paste("Processing", ontology_name, "...\n"))

    # Read GO enrichment results
    go_data <- read.csv(go_file)

    # Convert p-values to numeric and filter
    go_data$classicFisher <- as.numeric(go_data$classicFisher)
    go_data <- go_data %>%
        filter(!is.na(classicFisher), classicFisher < 0.05) %>%
        arrange(classicFisher) %>%
        head(10) %>%
        mutate(
            neg_log10_p = -log10(classicFisher),
            Term_short = ifelse(nchar(Term) > 50,
                               paste0(substr(Term, 1, 47), "..."),
                               Term),
            Term_with_count = paste0(Term_short, " (n=", Significant, ")"),
            Term_with_count = factor(Term_with_count, levels = rev(Term_with_count))
        )

    if (nrow(go_data) == 0) {
        cat(paste("  No significant terms found for", ontology_name, "\n\n"))
        return(NULL)
    }

    # Create bar chart - RED for L3
    p <- ggplot(go_data, aes(x = Term_with_count, y = neg_log10_p)) +
        geom_bar(stat = "identity", fill = "#D32F2F", alpha = 0.8) +
        coord_flip() +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
        labs(
            title = paste0("Top 10 Enriched GO Terms: ", ontology_name),
            subtitle = "L3-Upregulated Genes (N. americanus)",
            x = "",
            y = "-log10(p-value)"
        ) +
        theme_bw(base_size = 12) +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
            axis.text.y = element_text(size = 10)
        )

    # Save plot
    ggsave(paste0(output_name, ".pdf"), p, width = 10, height = 6, dpi = 300)
    ggsave(paste0(output_name, ".png"), p, width = 10, height = 6, dpi = 300)

    cat(paste("  ✓ Saved:", output_name, "\n\n"))

    return(go_data)
}

# Create plots for each ontology
bp_data <- plot_top_go_terms(
    "GO_enrichment_L3_upregulated_BP.csv",
    "Biological Process",
    "top10_GO_BP_L3_upregulated"
)

mf_data <- plot_top_go_terms(
    "GO_enrichment_L3_upregulated_MF.csv",
    "Molecular Function",
    "top10_GO_MF_L3_upregulated"
)

cc_data <- plot_top_go_terms(
    "GO_enrichment_L3_upregulated_CC.csv",
    "Cellular Component",
    "top10_GO_CC_L3_upregulated"
)

# Save summary tables
if (!is.null(bp_data)) {
    write.csv(bp_data[, c("GO.ID", "Term", "Significant", "Expected", "classicFisher")],
              "top10_GO_BP_L3_upregulated_table.csv", row.names = FALSE)
}

if (!is.null(mf_data)) {
    write.csv(mf_data[, c("GO.ID", "Term", "Significant", "Expected", "classicFisher")],
              "top10_GO_MF_L3_upregulated_table.csv", row.names = FALSE)
}

if (!is.null(cc_data)) {
    write.csv(cc_data[, c("GO.ID", "Term", "Significant", "Expected", "classicFisher")],
              "top10_GO_CC_L3_upregulated_table.csv", row.names = FALSE)
}

cat("================================================================\n")
cat("✓ All L3 GO term bar charts created!\n")
cat("  Files created:\n")
cat("    - top10_GO_BP_L3_upregulated.pdf/png\n")
cat("    - top10_GO_MF_L3_upregulated.pdf/png\n")
cat("    - top10_GO_CC_L3_upregulated.pdf/png\n")
cat("    - top10_GO_*_table.csv (data tables)\n")
cat("================================================================\n")
