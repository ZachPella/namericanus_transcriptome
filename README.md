# Necator americanus Transcriptome Analysis
## Differential Gene Expression Across Life Stages

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

---

## Overview

This repository contains the complete bioinformatics pipeline for analyzing differential gene expression in *Necator americanus* (human hookworm) across three life stages: L3 larvae, adult females, and adult males.

**Key Findings:**
- 15,990 genes quantified after quality filtering
- 9,116 differentially expressed genes between adult males and L3 larvae
- 8,376 differentially expressed genes between adult females and L3 larvae  
- 4,342 differentially expressed genes between adult females and males
<img width="652" height="968" alt="Untitled Diagram drawio (16)" src="https://github.com/user-attachments/assets/e5e89309-d680-4aaf-aadb-5e896a3edc79" />

---

## Table of Contents

- [Pipeline Overview](#pipeline-overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Directory Structure](#directory-structure)
- [Results](#results)
- [Citation](#citation)
- [Contact](#contact)

---

## Pipeline Overview
```
Raw FASTQ → QC → Alignment → Gene Prediction → Quantification → 
Differential Expression → Functional Annotation → GO Enrichment
```

**Analysis Steps:**
1. **Preprocessing**: Read QC, trimming, and alignment (HISAT2)
2. **Annotation**: Gene prediction using BRAKER3 with Augustus
3. **Quantification**: Read counting with featureCounts
4. **Differential Expression**: DESeq2 analysis (3 pairwise comparisons)
5. **Functional Annotation**: InterProScan domain/GO term annotation
6. **Enrichment Analysis**: Gene Ontology enrichment with topGO

---

## Requirements

### Software Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| fastp | 0.23.4 | Read trimming and QC |
| FastQC | 0.12.1 | Quality control reports |
| HISAT2 | 2.2.1 | RNA-seq alignment |
| SAMtools | 1.18 | BAM file processing |
| BRAKER3 | 3.0.8 | Gene prediction |
| featureCounts | 2.1.1 | Read quantification |
| R | 4.3.0+ | Statistical analysis |
| DESeq2 | 1.42.0 | Differential expression |
| topGO | 2.54.0 | GO enrichment |
| InterProScan | 5.69 | Functional annotation |

### R Packages
```r
install.packages(c("dplyr", "ggplot2", "pheatmap", "RColorBrewer", "tibble"))
BiocManager::install(c("DESeq2", "topGO"))
```

---

## Installation
```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/namericanus_transcriptome.git
cd namericanus_transcriptome

# Make scripts executable
chmod +x scripts/*.sh
```

---

## Usage

### Quick Start
```bash
# 1. Preprocessing (run in order)
sbatch scripts/preprocess_01_concatenate.sh
sbatch scripts/preprocess_02_qc_fastp.sh
sbatch scripts/preprocess_03_qc_fastqc.sh
sbatch scripts/preprocess_04_index_hisat2.sh
sbatch scripts/preprocess_05_align.sh
sbatch scripts/preprocess_06_bam.sh

# 2. Gene annotation
sbatch scripts/run_braker.sh

# 3. Differential expression analysis
sbatch scripts/1_feature_counts.sh
sbatch scripts/2_submit_deseq2.sh

# 4. Functional annotation
sbatch scripts/3_interproscan.sh
Rscript scripts/4_parse_IP.R
Rscript scripts/5_add_annotations_toDEseq.R

# 5. Enrichment and visualization
Rscript scripts/6_go_analysis.R
Rscript scripts/7_volcano_plots_all_comparisons.R
```

### Detailed Instructions

See [docs/PIPELINE.md](docs/PIPELINE.md) for step-by-step instructions.

---

## Directory Structure
```
.
├── scripts/                    # Analysis scripts
│   ├── preprocess_*.sh        # Read QC and alignment
│   ├── run_braker.sh          # Gene prediction
│   ├── 1_feature_counts.sh    # Quantification
│   ├── 2_submit_deseq2.sh     # DE analysis wrapper
│   ├── DESeq2_hookworm_analysis.R  # DESeq2 script
│   ├── 3_interproscan.sh      # Functional annotation
│   ├── 4_parse_IP.R           # Parse InterProScan
│   ├── 5_add_annotations_toDEseq.R  # Merge annotations
│   ├── 6_go_analysis.R        # GO enrichment
│   └── 7_volcano_plots_all_comparisons.R  # Visualization
├── docs/                      # Documentation
│   ├── PIPELINE.md           # Detailed pipeline guide
│   ├── METHODS.md            # Methods for publication
│   └── FLOWCHART.md          # Visual pipeline overview
├── results/                   # Example results (not included)
│   ├── figures/              # Publication figures
│   └── tables/               # Summary statistics
├── data/                     # Metadata
│   ├── sample_metadata.csv   # Sample information
│   └── README.md            # Data access information
├── README.md                 # This file
├── LICENSE                   # MIT License
└── environment.yml           # Conda environment (optional)
```

---

## Results

### Sample Information

| Life Stage | Replicates | Description |
|------------|------------|-------------|
| L3 larvae | 4 | Third-stage infective larvae |
| Adult Female | 4 | Mature female worms |
| Adult Male | 4 | Mature male worms |

**Sequencing:** Element Aviti platform, paired-end 150bp

### Key Statistics

- **Total reads processed:** ~XXX million
- **Average mapping rate:** XX%
- **Genes quantified:** 15,990
- **DE genes (padj < 0.05, |log2FC| > 1):**
  - Female vs L3: 8,376 genes
  - Male vs L3: 9,116 genes
  - Female vs Male: 4,342 genes

### Output Files

Generated files include:
- Differential expression results (CSV)
- Annotated gene lists with GO terms
- PCA plots
- Volcano plots
- Sample correlation heatmaps
- GO enrichment results

---

## Citation

If you use this pipeline or data, please cite:
```
[Author names]. (2025). Differential gene expression analysis across life stages 
of Necator americanus. [Journal name]. DOI: XX.XXXX/XXXXXX
```

**Raw data availability:** [SRA/GEO accession number]

---



## Contact

**Zachary Pella**  
APHL-CDC Bioinformatics Fellow  
University of Nebraska Medical Center  
Nebraska Public Health Laboratory  
Email: [your.email@unmc.edu]

**PI:** [Dr. Name]  
Email: [pi.email@unmc.edu]

---

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- APHL-CDC Fellowship Program
- Dr. [Your PI's name] and the Fauver Lab
- UNMC Nebraska Public Health Laboratory
- UNMC High Performance Computing Core

---

**Last Updated:** December 2024
