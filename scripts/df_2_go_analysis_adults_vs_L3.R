#!/usr/bin/env Rscript
# ==============================================================================
# GO Enrichment Analysis: Adult vs L3
# ==============================================================================
library(topGO)
library(dplyr)
library(tibble)
BASE_DIR <- "/work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data"
OUTPUT_DIR <- paste0(BASE_DIR, "/differential_expression_total_adults_vs_L3")
setwd(OUTPUT_DIR)
cat("================================================================\n")
cat("GO Enrichment Analysis: Adult vs L3\n")
cat("================================================================\n\n")
# Load annotations
annotations <- read.csv(paste0(BASE_DIR, "/functional_annotation/gene_annotations_master.csv"))
# Clean GO terms - remove (InterPro) suffix
annotations$GO_terms <- gsub("\\(InterPro\\)", "", annotations$GO_terms)
cat(paste("Total genes in annotation database:", nrow(annotations), "\n"))
cat(paste("Genes with GO terms:", sum(!is.na(annotations$GO_terms)), "\n\n"))
# Create gene2GO mapping for topGO
gene2GO <- annotations %>%
    filter(!is.na(GO_terms)) %>%
    select(gene_id, GO_terms) %>%
    mutate(GO_list = strsplit(GO_terms, ";")) %>%
    select(gene_id, GO_list) %>%
    deframe()
# Load significant genes
sig_genes <- read.csv("DESeq2_significant_Adult_vs_L3.csv")
cat(paste("Significant DEGs:", nrow(sig_genes), "\n"))
cat(paste("  Upregulated in Adult:", sum(sig_genes$log2FoldChange > 0), "\n"))
cat(paste("  Downregulated in Adult (= L3 upregulated):", sum(sig_genes$log2FoldChange < 0), "\n\n"))
# Function to run GO enrichment
run_go_enrichment <- function(sig_genes, all_genes, comp_name) {
    cat(paste("Running GO enrichment for", comp_name, "...\n\n"))
    # Create gene list (1 = significant, 0 = not significant)
    geneList <- factor(as.integer(all_genes %in% sig_genes))
    names(geneList) <- all_genes
    # Build topGO object for each ontology
    for (ont in c("BP", "MF", "CC")) {
        cat(paste("  Ontology:", ont, "\n"))
        GOdata <- new("topGOdata",
                      ontology = ont,
                      allGenes = geneList,
                      annot = annFUN.gene2GO,
                      gene2GO = gene2GO)
        # Run Fisher's exact test
        resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
        # Get results table
        allRes <- GenTable(GOdata,
                          classicFisher = resultFisher,
                          orderBy = "classicFisher",
                          ranksOf = "classicFisher",
                          topNodes = 50)
        # Save results
        write.csv(allRes,
                  paste0("GO_enrichment_", comp_name, "_", ont, ".csv"),
                  row.names = FALSE)
        cat(paste("    Significant terms (p<0.05):",
                  sum(as.numeric(allRes$classicFisher) < 0.05), "\n"))
    }
    cat("\n")
}

# ==============================================================================
# 1. ALL significant DEGs (original analysis)
# ==============================================================================
all_genes <- annotations$gene_id
run_go_enrichment(sig_genes$gene_id, all_genes, "Adult_vs_L3")

# ==============================================================================
# 2. ADULT-UPREGULATED genes only (log2FC > 0)
# ==============================================================================
cat("================================================================\n")
cat("GO Enrichment: ADULT-UPREGULATED genes\n")
cat("================================================================\n\n")

adult_upregulated <- sig_genes %>%
    filter(log2FoldChange > 0) %>%
    pull(gene_id)

cat(paste("Adult-upregulated genes:", length(adult_upregulated), "\n\n"))
run_go_enrichment(adult_upregulated, all_genes, "Adult_upregulated")

# ==============================================================================
# 3. L3-UPREGULATED genes only (log2FC < 0)
# ==============================================================================
cat("================================================================\n")
cat("GO Enrichment: L3-UPREGULATED genes\n")
cat("================================================================\n\n")

l3_upregulated <- sig_genes %>%
    filter(log2FoldChange < 0) %>%
    pull(gene_id)

cat(paste("L3-upregulated genes:", length(l3_upregulated), "\n\n"))
run_go_enrichment(l3_upregulated, all_genes, "L3_upregulated")

cat("================================================================\n")
cat("✓ GO enrichment analysis complete!\n")
cat("  Files created:\n")
cat("    - GO_enrichment_Adult_vs_L3_*.csv (all DEGs)\n")
cat("    - GO_enrichment_Adult_upregulated_*.csv (Adult up)\n")
cat("    - GO_enrichment_L3_upregulated_*.csv (L3 up)\n")
cat("================================================================\n")
