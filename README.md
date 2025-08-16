# bulk-rnaseq-deseq2

Reproducible **bulk RNA-seq differential expression** example using **DESeq2** on a small synthetic dataset (10 genes × 12 samples). 
This is a clean, minimal repo you can show on your CV/GitHub to demonstrate RNA-seq analysis skills.

## What this repo shows
- Loading counts + metadata
- Running DESeq2 differential expression (infected vs control)
- Volcano plot + MA plot
- Export of results table as CSV

## Quick start

### 1) Setup R environment
You need R (≥4.2) and these packages:
```
install.packages(c("BiocManager"))
BiocManager::install(c("DESeq2", "apeglm"))
install.packages(c("ggplot2", "readr", "dplyr"))
```

### 2) Run the pipeline
From the repo root:
```
Rscript scripts/deseq2_analysis.R   --counts data/counts.csv   --metadata data/metadata.csv   --condition_col condition   --contrast infected control   --outdir results
```

### 3) Outputs
- `results/deseq2_results.csv` — full results table (log2FC, pvalue, padj)
- `results/volcano_plot.png`
- `results/ma_plot.png`
- `results/sessionInfo.txt`

## Data
Synthetic data in `data/`:
- `counts.csv` — gene counts (rows=genes, cols=samples)
- `metadata.csv` — sample metadata with `sample` and `condition`

## Notes
This is a demo with tiny data to keep things light. On real projects you'd add:
- QC (FastQC/MultiQC), alignment (STAR/HISAT2), featureCounts/Salmon outputs
- Snakemake/Nextflow workflow
- Better covariate handling (batch, sex, etc.)
