#!/usr/bin/env Rscript

library(DESeq2)
library(ggplot2)
library(apeglm)

args <- commandArgs(trailingOnly = TRUE)

counts_data <- read.csv("countmatrix.csv", header=TRUE, row.names="X")

colData <- read.csv(args[1], header=TRUE, row.names="sample")

tryCatch({

  dds <- DESeqDataSetFromMatrix(countData = counts_data,
                                colData = colData,
                                design = ~ condition)

  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]


  dds$condition <- relevel(dds$condition, ref="WT")

  dds <- DESeq(dds)

  res <- results(dds, alpha = 0.05)

  png("deseq2_output/MA_plot.png", width = 800, height = 600)
  plotMA(res)
  dev.off()

  vsd <- vst(dds, blind = TRUE)
  png("deseq2_output/PCA_plot.png", width = 800, height = 600)
  plotPCA(vsd, intgroup = "condition")
  dev.off()

  resLFC <- lfcShrink(dds, coef = "condition_snf2_vs_WT", type = "apeglm", res = res)

  png("deseq2_output/MA_shrunkenLFC_plot.png", width = 800, height = 600)
  plotMA(resLFC, ylim = c(-8, 8))
  dev.off()

  png("deseq2_output/compare_shrunkenLFC_plot.png", width = 800, height = 600)
  plot(
    res$log2FoldChange,
    resLFC$log2FoldChange,
    xlab = "Unshrunken LFC",
    ylab = "Shrunken LFC"
  )
  abline(0, 1, col = "red")
  dev.off()

  resdf <- as.data.frame(resLFC)
  resdf$sig <- !is.na(resdf$padj) & 
              resdf$padj < 0.05 & 
              abs(resdf$log2FoldChange) > 1

  volcano <- ggplot(resdf, aes(log2FoldChange, -log10(pvalue), colour = sig)) +
    geom_point(alpha = 0.4, size = 0.8) +
    scale_colour_manual(values = c("grey70", "#c0392b")) +
    theme_minimal() + labs(colour = "padj<0.05 & |LFC|>1")
  ggsave("deseq2_output/Volcano_plot.png", plot = volcano)

  write.csv(as.data.frame(resLFC[order(resLFC$padj), ]),
            "deseq2_output/de_snf2_vs_WT.csv") },
error = function(e) {
  cat("An error occured: ", conditionMessage(e), "\n")
  NA
}, finally = {
  cat("DeSEQ2 has been executed.\n")
})
