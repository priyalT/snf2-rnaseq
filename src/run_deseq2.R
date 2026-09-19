library(DESeq2)
library(ggplot2)
library(apeglm)
library(tidyverse)

counts_data <- read.csv("../results/countmatrix.csv", header=TRUE, row.names="X")
head(counts_data)

colData <- read.csv("../data/samplesheet.csv", header=TRUE, row.names="sample")

all(colnames(counts_data) %in% rownames(colData))
all(colnames(counts_data) == rownames(colData))

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ condition)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds

dds$condition <- relevel(dds$condition, ref="WT")

dds$condition    

dds <- DESeq(dds)

res <- results(dds, alpha = 0.05)
summary(res)

png("../plots/MA_plot.png", width = 800, height = 600)
plotMA(res)
dev.off()

vsd <- vst(dds, blind = TRUE)
png("../plots/PCA_plot.png", width = 800, height = 600)
plotPCA(vsd, intgroup = "condition")
dev.off()

resultsNames(dds)
resLFC <- lfcShrink(dds, coef = "condition_snf2_vs_WT", type = "apeglm", res = res)

png("../plots/MA_shrunkenLFC_plot.png", width = 800, height = 600)
plotMA(resLFC)
dev.off()

png("../plots/compare_shrunkenLFC_plot.png", width = 800, height = 600)
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
ggsave("../plots/Volcano_plot.png", plot = volcano)

write.csv(as.data.frame(resLFC[order(resLFC$padj), ]),
          "../results/de_snf2_vs_WT.csv")