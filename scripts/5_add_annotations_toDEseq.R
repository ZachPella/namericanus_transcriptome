#!/usr/bin/env Rscript

# ==============================================================================
# Add Functional Annotations to DESeq2 Results
# ==============================================================================

library(dplyr)

BASE_DIR <- "/work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data"
ANNOT_FILE <- paste0(BASE_DIR, "/functional_annotation/gene_annotations_master.csv")
DE_DIR <- paste0(BASE_DIR, "/differential_expression")

setwd(DE_DIR)

cat("================================================================\n")
cat("Adding Functional Annotations to DESeq2 Results\n")
cat("================================================================\n\n")

# Load annotations
annotations <- read.csv(ANNOT_FILE)

cat(paste("Loaded annotations for", nrow(annotations), "genes\n\n"))

# Process each comparison
comparisons <- c("Female_vs_L3", "Male_vs_L3", "Female_vs_Male")

for (comp in comparisons) {
    cat(paste("Processing", comp, "...\n"))
    
    # Read DESeq2 results
    results_file <- paste0("DESeq2_results_", comp, ".csv")
    sig_file <- paste0("DESeq2_significant_", comp, ".csv")
    
    if (!file.exists(results_file)) {
        cat(paste("  Warning: File not found:", results_file, "\n"))
        next
    }
    
    # Read and merge
    results <- read.csv(results_file)
    results_annotated <- results %>%
        left_join(annotations, by = "gene_id")
    
    # Save annotated results
    write.csv(results_annotated,
              paste0("DESeq2_results_ANNOTATED_", comp, ".csv"),
              row.names = FALSE)
    
    # Also annotate significant genes
    if (file.exists(sig_file)) {
        sig <- read.csv(sig_file)
        sig_annotated <- sig %>%
            left_join(annotations, by = "gene_id")
        
        write.csv(sig_annotated,
                  paste0("DESeq2_significant_ANNOTATED_", comp, ".csv"),
                  row.names = FALSE)
        
        # Summary
        n_sig <- nrow(sig_annotated)
        n_with_go <- sum(!is.na(sig_annotated$GO_terms))
        n_with_domain <- sum(!is.na(sig_annotated$InterPro_IDs))
        
        cat(paste("  Significant genes:", n_sig, "\n"))
        cat(paste("    With GO terms:", n_with_go, 
                  "(", round(100*n_with_go/n_sig, 1), "%)\n"))
        cat(paste("    With domains:", n_with_domain,
                  "(", round(100*n_with_domain/n_sig, 1), "%)\n\n"))
    }
}

cat("================================================================\n")
cat("✓ Annotation complete!\n")
cat("================================================================\n")
