# Detailed Pipeline Guide

## Step-by-Step Instructions

### Prerequisites
- SLURM job scheduler
- All software installed (see README.md)
- Sufficient storage space (~500GB)
- At least 88GB RAM for InterProScan

---

## Phase 1: Preprocessing

### Step 1: Concatenate Read Files
```bash
sbatch scripts/preprocess_01_concatenate.sh
```
**Purpose:** Combine technical replicates  
**Time:** ~30 minutes  
**Output:** `concatenated_reads/`

### Step 2: Quality Trimming
```bash
sbatch scripts/preprocess_02_qc_fastp.sh
```
**Purpose:** Adapter trimming, quality filtering  
**Time:** ~2 hours  
**Output:** `trimmed_reads/`

### Step 3: Quality Control Reports
```bash
sbatch scripts/preprocess_03_qc_fastqc.sh
```
**Purpose:** Generate QC metrics  
**Time:** ~1 hour  
**Output:** `fastqc_reports/`

### Step 4: Index Reference Genome
```bash
sbatch scripts/preprocess_04_index_hisat2.sh
```
**Purpose:** Build HISAT2 index  
**Time:** ~30 minutes  
**Output:** `reference/*.ht2`

### Step 5: Alignment
```bash
sbatch scripts/preprocess_05_align.sh
```
**Purpose:** Align reads to reference  
**Time:** ~4 hours  
**Output:** `sam_files/`

### Step 6: BAM Processing
```bash
sbatch scripts/preprocess_06_bam.sh
```
**Purpose:** Convert to BAM, sort, index  
**Time:** ~2 hours  
**Output:** `bam_files/*.sorted.bam`

---

## Phase 2: Gene Annotation

### Step 7: BRAKER3 Gene Prediction
```bash
sbatch scripts/run_braker.sh
```
**Purpose:** Predict gene structures  
**Time:** ~24-48 hours  
**Output:** `braker_annotation/braker_output/Augustus/augustus.hints.gtf`

**Check output:**
```bash
grep -c "gene_id" braker_annotation/braker_output/Augustus/augustus.hints.gtf
```

---

## Phase 3: Quantification & Differential Expression

### Step 8: featureCounts
```bash
sbatch scripts/1_feature_counts.sh
```
**Purpose:** Count reads per gene  
**Time:** ~1 hour  
**Output:** `differential_expression/gene_counts_all_samples.txt`

**Verify counts:**
```bash
head differential_expression/gene_counts_all_samples.txt.summary
```

### Step 9: DESeq2 Analysis
```bash
sbatch scripts/2_submit_deseq2.sh
```
**Purpose:** Differential expression testing  
**Time:** ~30 minutes  
**Output:**
- `DESeq2_results_*.csv`
- `DESeq2_significant_*.csv`
- `PCA_plot.pdf`
- `volcano_*.pdf`

**Check results:**
```bash
wc -l differential_expression/DESeq2_significant_*.csv
```

---

## Phase 4: Functional Annotation

### Step 10: InterProScan
```bash
sbatch scripts/3_interproscan.sh
```
**Purpose:** Annotate protein domains and functions  
**Time:** ~24-48 hours  
**Memory:** 88GB required  
**Output:** `functional_annotation/interproscan/*.tsv`

### Step 11: Parse InterProScan Results
```bash
Rscript scripts/4_parse_IP.R
```
**Purpose:** Convert to gene-level annotations  
**Time:** ~5 minutes  
**Output:** `functional_annotation/gene_annotations_master.csv`

### Step 12: Merge Annotations with DESeq2
```bash
Rscript scripts/5_add_annotations_toDEseq.R
```
**Purpose:** Add GO terms to DE results  
**Time:** ~2 minutes  
**Output:** `differential_expression/DESeq2_*_ANNOTATED_*.csv`

---

## Phase 5: Enrichment & Visualization

### Step 13: GO Enrichment Analysis
```bash
Rscript scripts/6_go_analysis.R
```
**Purpose:** Identify enriched biological processes  
**Time:** ~10 minutes  
**Output:** `differential_expression/GO_enrichment_*.csv`

### Step 14: Enhanced Volcano Plots
```bash
Rscript scripts/7_volcano_plots_all_comparisons.R
```
**Purpose:** Publication-quality figures  
**Time:** ~2 minutes  
**Output:** `differential_expression/volcano_plot_*.pdf`

---

## Verification Checklist

After completing the pipeline, verify:

- [ ] All 12 BAM files generated and indexed
- [ ] BRAKER GTF contains gene predictions
- [ ] featureCounts summary shows >70% assignment rate
- [ ] PCA plot shows clear separation by life stage
- [ ] Significant genes identified in each comparison
- [ ] GO terms assigned to >50% of genes
- [ ] Volcano plots generated for all comparisons

---

## Troubleshooting

### Common Issues

**Problem:** featureCounts fails with "GTF file not found"  
**Solution:** Verify BRAKER completed successfully and check path in script

**Problem:** DESeq2 gives "samples and colData different"  
**Solution:** Check sample names match between counts and metadata

**Problem:** InterProScan runs out of memory  
**Solution:** Increase to 88GB RAM or split protein file into batches

**Problem:** GO enrichment finds no significant terms  
**Solution:** Check that GO terms were successfully parsed from InterProScan

---

## Performance Notes

**Total pipeline time:** ~4-5 days  
**Storage required:** ~500GB  
**Peak RAM usage:** 88GB (InterProScan)  
**Recommended compute:** 16-32 cores for parallel steps

