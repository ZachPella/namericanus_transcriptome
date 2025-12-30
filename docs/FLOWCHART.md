# Pipeline Flowchart

See the ASCII flowchart in the main README or view the visual version in `docs/figures/pipeline_flowchart.png` (if available).

## Pipeline Stages

1. **Preprocessing** (Steps 1-6)
2. **Annotation** (Step 7)
3. **Quantification** (Steps 8-9)
4. **Functional Annotation** (Steps 10-12)
5. **Enrichment** (Steps 13-14)

## Data Flow
```
Raw FASTQ Files (12 samples)
    ↓
Trimmed & QC'd Reads
    ↓
Aligned Reads (BAM files)
    ↓
Gene Predictions (BRAKER3 GTF)
    ↓
Gene Count Matrix (15,990 genes × 12 samples)
    ↓
Differential Expression Results (3 comparisons)
    ↓
Functional Annotations (GO terms, domains)
    ↓
Enriched Biological Processes
```
