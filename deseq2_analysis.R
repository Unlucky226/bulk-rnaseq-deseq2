#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(DESeq2)
  library(ggplot2)
})

option_list <- list(
  make_option(c("--counts"), type="character", help="Path to counts.csv"),
  make_option(c("--metadata"), type="character", help="Path to metadata.csv"),
  make_option(c("--condition_col"), type="character", default="condition", help="Condition column name"),
  make_option(c("--contrast"), type="character", default="infected control", help="'case control' (space separated)"),
  make_option(c("--outdir"), type="character", default="results", help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

counts <- read.csv(opt$counts, row.names=1, check.names=FALSE)
meta <- read.csv(opt$metadata, stringsAsFactors=FALSE)
stopifnot("sample" %in% colnames(meta))
rownames(meta) <- meta$sample
meta[[opt$condition_col]] <- factor(meta[[opt$condition_col]])

# Ensure columns match
common <- intersect(colnames(counts), rownames(meta))
counts <- counts[, common, drop=FALSE]
meta <- meta[common, , drop=FALSE]

dds <- DESeqDataSetFromMatrix(countData=round(as.matrix(counts)),
                              colData=meta,
                              design=as.formula(paste0("~ ", opt$condition_col)))
dds <- DESeq(dds)

parts <- strsplit(opt$contrast, " +")[[1]]
case <- parts[1]; ctrl <- parts[2]

res <- results(dds, contrast=c(opt$condition_col, case, ctrl))
res <- lfcShrink(dds, coef=2, type="apeglm", res=res)
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df %>% arrange(padj)
write.csv(res_df, file=file.path(opt$outdir, "deseq2_results.csv"), row.names=FALSE)

# MA plot
png(file.path(opt$outdir, "ma_plot.png"), width=1200, height=900, res=150)
plotMA(res, ylim=c(-4,4))
dev.off()

# Volcano plot
res_df$neglog10padj <- -log10(res_df$padj)
res_df$neglog10padj[!is.finite(res_df$neglog10padj)] <- NA

png(file.path(opt$outdir, "volcano_plot.png"), width=1200, height=900, res=150)
ggplot(res_df, aes(x=log2FoldChange, y=neglog10padj)) +
  geom_point() +
  geom_vline(xintercept=c(-1,1), linetype="dashed") +
  geom_hline(yintercept=1.3, linetype="dashed") +
  xlab("log2 Fold Change") + ylab("-log10 adjusted p-value") +
  ggtitle("Volcano plot")
dev.off()

# Session info
sink(file.path(opt$outdir, "sessionInfo.txt"))
print(sessionInfo())
sink()
