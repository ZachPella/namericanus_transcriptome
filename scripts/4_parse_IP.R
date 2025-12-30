#!/usr/bin/env Rscript

# ==============================================================================
# Parse InterProScan Results and Create Gene Annotation Database
# ==============================================================================

library(dplyr)
library(tidyr)

BASE_DIR <- "/work/fauverlab/zachpella/hookworm/transcriptome_Dec2025_raw_aviti_data"
INTERPRO_FILE <- paste0(BASE_DIR, "/functional_annotation/interproscan/augustus.hints.aa.tsv")
OUTPUT_DIR <- paste0(BASE_DIR, "/functional_annotation")

cat("================================================================\n")
cat("Parsing InterProScan Results\n")
cat("================================================================\n\n")

# Read InterProScan TSV
# Columns: 1=protein_id, 2=md5, 3=length, 4=analysis, 5=signature, 6=description, 
#          7=start, 8=end, 9=score, 10=status, 11=date, 12=(InterPro ID), 
#          13=(InterPro description), 14=GO terms, 15=Pathways

interpro <- read.delim(INTERPRO_FILE, header = FALSE, stringsAsFactors = FALSE)
colnames(interpro) <- c("protein_id", "md5", "length", "analysis", "signature", 
                        "description", "start", "end", "score", "status", "date",
                        "interpro_id", "interpro_desc", "go_terms", "pathways")

cat(paste("Total annotations:", nrow(interpro), "\n"))
cat(paste("Unique proteins:", length(unique(interpro$protein_id)), "\n\n"))

# ==============================================================================
# 1. Extract Gene ID from protein ID
# ==============================================================================
# BRAKER protein IDs are like: g1234.t1, g1234.t2, etc.
# We need to extract the gene ID (g1234)

interpro <- interpro %>%
    mutate(gene_id = sub("\\.t.*$", "", protein_id))

# ==============================================================================
# 2. Create Gene-Level GO Term Mapping
# ==============================================================================

cat("Creating gene-level GO term database...\n")

go_mapping <- interpro %>%
    filter(!is.na(go_terms), go_terms != "-") %>%
    select(gene_id, go_terms) %>%
    separate_rows(go_terms, sep = "\\|") %>%
    distinct() %>%
    group_by(gene_id) %>%
    summarize(
        GO_terms = paste(sort(unique(go_terms)), collapse = ";"),
        GO_count = n()
    )

cat(paste("Genes with GO terms:", nrow(go_mapping), "\n\n"))

# ==============================================================================
# 3. Create Gene-Level Domain/Family Mapping
# ==============================================================================

cat("Creating domain annotation database...\n")

domain_mapping <- interpro %>%
    filter(!is.na(interpro_id), interpro_id != "-") %>%
    select(gene_id, interpro_id, interpro_desc) %>%
    distinct() %>%
    group_by(gene_id) %>%
    summarize(
        InterPro_IDs = paste(unique(interpro_id), collapse = ";"),
        InterPro_descriptions = paste(unique(interpro_desc), collapse = ";"),
        Domain_count = n()
    )

cat(paste("Genes with InterPro domains:", nrow(domain_mapping), "\n\n"))

# ==============================================================================
# 4. Create Gene-Level Pathway Mapping
# ==============================================================================

cat("Creating pathway annotation database...\n")

pathway_mapping <- interpro %>%
    filter(!is.na(pathways), pathways != "-") %>%
    select(gene_id, pathways) %>%
    separate_rows(pathways, sep = "\\|") %>%
    distinct() %>%
    group_by(gene_id) %>%
    summarize(
        Pathways = paste(sort(unique(pathways)), collapse = ";"),
        Pathway_count = n()
    )

cat(paste("Genes with pathway annotations:", nrow(pathway_mapping), "\n\n"))

# ==============================================================================
# 5. Create Master Annotation Table
# ==============================================================================

cat("Creating master gene annotation table...\n")

# Get all unique genes
all_genes <- data.frame(gene_id = unique(interpro$gene_id))

# Merge all annotations
gene_annotations <- all_genes %>%
    left_join(go_mapping, by = "gene_id") %>%
    left_join(domain_mapping, by = "gene_id") %>%
    left_join(pathway_mapping, by = "gene_id") %>%
    mutate(
        has_GO = !is.na(GO_terms),
        has_domain = !is.na(InterPro_IDs),
        has_pathway = !is.na(Pathways)
    )

# Save master annotation table
write.csv(gene_annotations, 
          paste0(OUTPUT_DIR, "/gene_annotations_master.csv"),
          row.names = FALSE)

cat(paste("Master annotation table saved:", nrow(gene_annotations), "genes\n\n"))

# ==============================================================================
# 6. Summary Statistics
# ==============================================================================

cat("================================================================\n")
cat("ANNOTATION SUMMARY\n")
cat("================================================================\n\n")

cat(paste("Total genes annotated:", nrow(gene_annotations), "\n"))
cat(paste("Genes with GO terms:", sum(gene_annotations$has_GO), 
          "(", round(100*sum(gene_annotations$has_GO)/nrow(gene_annotations), 1), "%)\n"))
cat(paste("Genes with InterPro domains:", sum(gene_annotations$has_domain),
          "(", round(100*sum(gene_annotations$has_domain)/nrow(gene_annotations), 1), "%)\n"))
cat(paste("Genes with pathway annotations:", sum(gene_annotations$has_pathway),
          "(", round(100*sum(gene_annotations$has_pathway)/nrow(gene_annotations), 1), "%)\n\n"))

cat("================================================================\n")
