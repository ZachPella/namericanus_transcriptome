# Methods Section for Publication

## RNA Sequencing and Library Preparation

Total RNA was extracted from *Necator americanus* samples representing three life stages: third-stage infective larvae (L3, n=4), adult females (n=4), and adult males (n=4). RNA quality was assessed using [method]. Libraries were prepared using [kit name] and sequenced on the Element Aviti platform, generating paired-end 150bp reads.

## Read Processing and Quality Control

Raw sequencing reads were processed using fastp v0.23.4 for adapter trimming and quality filtering (minimum quality score: 20, minimum length: 50bp). Quality metrics were assessed using FastQC v0.12.1. Processed reads were aligned to the *N. americanus* reference genome (GenBank assembly: [accession]) using HISAT2 v2.2.1 with default parameters for paired-end RNA-seq data.

## Gene Annotation

Gene structures were predicted using BRAKER3 v3.0.8 in RNA-seq mode with Augustus as the gene predictor. The pipeline integrated RNA-seq alignment evidence to generate gene models. The resulting annotation contained [N] protein-coding genes.

## Differential Expression Analysis

Gene-level read counts were quantified using featureCounts v2.1.1 (Subread package) with the following parameters: paired-end mode (-p), fragment counting (-B), exclusion of chimeric fragments (-C), and primary alignments only (--primary). 

Differential expression analysis was performed using DESeq2 v1.42.0 in R v4.3.0. Genes with fewer than 10 total counts across all samples were filtered prior to analysis. The experimental design included life stage as the primary factor with three levels (L3, Female, Male). Size factors were calculated using the median-of-ratios method, and dispersion estimates were obtained using the default parametric fit.

Three pairwise comparisons were performed:
1. Adult females vs. L3 larvae
2. Adult males vs. L3 larvae  
3. Adult females vs. Adult males

Differentially expressed genes (DEGs) were identified using an adjusted p-value threshold of 0.05 (Benjamini-Hochberg FDR correction) and an absolute log2 fold change > 1 (≥2-fold change). Variance-stabilizing transformation (VST) was applied for visualization purposes including principal component analysis and hierarchical clustering.

## Functional Annotation

Predicted protein sequences were functionally annotated using InterProScan v5.69, querying the following databases: Pfam, SMART, SUPERFAMILY, Gene3D, PRINTS, ProSiteProfiles, and HAMAP. Gene Ontology (GO) terms and pathway annotations were extracted and mapped to gene identifiers.

## Gene Ontology Enrichment Analysis

GO enrichment analysis was performed using the topGO package v2.54.0 in R. Over-representation of GO terms in the set of differentially expressed genes was assessed using Fisher's exact test with the "classic" algorithm for each of the three GO ontologies (Biological Process, Molecular Function, Cellular Component). Terms with p-values < 0.05 were considered significantly enriched.

## Data Visualization

All visualizations were generated in R using ggplot2 v3.4.0 and pheatmap v1.0.12. Volcano plots display log2 fold changes versus -log10 adjusted p-values, with genes colored by significance category. Principal component analysis plots show sample relationships based on the top 500 most variable genes.

## Statistical Analysis

All statistical analyses were performed in R v4.3.0. Multiple testing correction was applied using the Benjamini-Hochberg method to control the false discovery rate. Results were considered statistically significant at FDR < 0.05 unless otherwise specified.

## Data Availability

Raw sequencing data have been deposited in the NCBI Sequence Read Archive (SRA) under BioProject accession [PRJNAXXXXXX]. Processed data, including gene count matrices and differential expression results, are available at [GitHub repository URL]. Analysis scripts are available at [GitHub repository URL].
