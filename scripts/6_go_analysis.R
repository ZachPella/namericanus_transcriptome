#!/usr/bin/env Rscript

# ==============================================================================
# GO Enrichment Analysis using topGO
# ==============================================================================

library(topGO)
library(dplyr)
library(tibble)

BASE_DIR <- "/work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data"
setwd(paste0(BASE_DIR, "/differential_expression"))

# Load annotations
annotations <- read.csv(paste0(BASE_DIR, "/functional_annotation/gene_annotations_master.csv"))

# Clean GO terms - remove (InterPro) suffix
annotations$GO_terms <- gsub("\\(InterPro\\)", "", annotations$GO_terms)

# Create gene2GO mapping for topGO
gene2GO <- annotations %>%
    filter(!is.na(GO_terms)) %>%
    select(gene_id, GO_terms) %>%
    mutate(GO_list = strsplit(GO_terms, ";")) %>%
    select(gene_id, GO_list) %>%
    deframe()

# Function to run GO enrichment
run_go_enrichment <- function(sig_genes, all_genes, comp_name) {
    cat(paste("Running GO enrichment for", comp_name, "...\n"))

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
}

# Run for each comparison
comparisons <- c("Female_vs_L3", "Male_vs_L3", "Female_vs_Male")

for (comp in comparisons) {
    sig_file <- paste0("DESeq2_significant_", comp, ".csv")

    if (file.exists(sig_file)) {
        sig_genes <- read.csv(sig_file)$gene_id
        all_genes <- annotations$gene_id

        run_go_enrichment(sig_genes, all_genes, comp)
    }
}

cat("\n✓ GO enrichment analysis complete!\n")
