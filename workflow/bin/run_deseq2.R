#!/usr/bin/env Rscript

library(DESeq2)
library(ggplot2)
library(apeglm)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep("^--file=", args)
  if (length(match) > 0) return(dirname(normalizePath(sub("^--file=", "", args[match]))))
  return(getwd())
}
theme_file <- file.path(get_script_dir(), "plot_theme.R")
if (file.exists(theme_file)) source(theme_file)

theme_set(custom_theme)

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
  resdf <- as.data.frame(res)
  resdf$gene <- rownames(resdf)
  resdf$sig <- !is.na(resdf$padj) & resdf$padj < 0.05
  
  top_genes <- head(resdf[order(resdf$padj), ], 10)

  ma_plot <- ggplot(resdf, aes(x = baseMean, y = log2FoldChange, colour = sig)) +
    geom_point(alpha = 0.4, size = 0.8) +
    scale_x_log10() +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    scale_colour_manual(values = c("grey70", "#0073C2FF"), labels = c("Not significant", "Significant")) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    geom_text(data = top_genes, aes(label = gene), size = 3, show.legend = FALSE, vjust = -0.5, check_overlap = TRUE) +
    labs(title = "MA Plot (Unshrunken)", x = "Mean of Normalized Counts", y = "Log2 Fold Change", colour = "Status")
  
  ggsave("deseq2_output/MA_plot.png", plot = ma_plot, width = 8, height = 6)
  
  vsd <- vst(dds, blind = FALSE)
  
  pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  
  pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
    geom_point(size = 3) +
    geom_text(size = 3, show.legend = FALSE, check_overlap = TRUE, vjust = -0.5) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    labs(
      title = "PCA: snf2 (KO) vs WT",
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance")
    )
  ggsave("deseq2_output/PCA_plot.png", plot = pca_plot, width = 8, height = 6)
  
  resLFC <- lfcShrink(dds, coef = "condition_snf2_vs_WT", type = "apeglm", res = res)
  
  shrunkenresdf <- as.data.frame(resLFC)
  shrunkenresdf$gene <- rownames(shrunkenresdf)
  shrunkenresdf$sig <- !is.na(shrunkenresdf$padj) & shrunkenresdf$padj < 0.05
  top_genes_shrunken <- head(shrunkenresdf[order(shrunkenresdf$padj), ], 10)
  shrunken_ma_plot <- ggplot(shrunkenresdf, aes(x = baseMean, y = log2FoldChange, colour = sig)) +
    geom_point(alpha = 0.4, size = 0.8) +
    scale_x_log10() +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    scale_colour_manual(values = c("grey70", "#0073C2FF"), labels = c("Not significant", "Significant")) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    geom_text(data = top_genes_shrunken, aes(label = gene), size = 3, show.legend = FALSE, vjust = -0.5, check_overlap = TRUE) +
    labs(title = "MA Plot (Shrunken)", x = "Mean of Normalized Counts", y = "Log2 Fold Change", colour = "Status")
  
  ggsave("deseq2_output/MA_shrunkenLFC_plot.png", plot = shrunken_ma_plot, width = 8, height = 6)
  
  
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

  resdf$gene <- rownames(resdf)
  top_genes_volcano <- head(resdf[order(resdf$padj), ], 10)
  volcano <- ggplot(resdf, aes(x = log2FoldChange, y = -log10(pvalue), colour = sig)) +
    geom_point(alpha = 0.4, size = 0.8) +
    scale_colour_manual(values = c("grey70", "#c0392b"), labels = c("Not significant", "Significant")) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.3) +
    geom_text(data = top_genes_volcano, aes(label = gene), size = 3, show.legend = FALSE, vjust = -0.5, check_overlap = TRUE) +
    labs(
      title = "Volcano Plot: snf2 (KO) vs WT",
      subtitle = "1,372 significant genes (padj < 0.05)",
      colour = "padj < 0.05 & |LFC| > 1"
    )
  ggsave("deseq2_output/Volcano_plot.png", plot = volcano, width = 8, height = 6)

  # Heatmap of top 50 differentially expressed genes (base R stats::heatmap)
  top50_genes <- head(rownames(resdf[order(resdf$padj), ]), 50)
  mat <- assay(vsd)[top50_genes, ]
  mat <- mat - rowMeans(mat)
  
  col_side_colors <- ifelse(colData$condition == "WT", "#4DAF4A", "#E41A1C")
  
  png("deseq2_output/Heatmap_top50.png", width = 800, height = 1000)
  heatmap(mat,
          scale = "none",
          ColSideColors = col_side_colors,
          col = colorRampPalette(c("#0073C2FF", "white", "#E41A1C"))(50),
          main = "Top 50 DE Genes (snf2 knockout vs WT)",
          margins = c(8, 8))
  dev.off()

  write.csv(as.data.frame(resLFC[order(resLFC$padj), ]),
            "deseq2_output/de_snf2_vs_WT.csv") },
error = function(e) {
  cat("An error occured: ", conditionMessage(e), "\n")
  quit(save = "no", status = 1)
}, finally = {
  cat("DeSEQ2 has been executed.\n")
})
