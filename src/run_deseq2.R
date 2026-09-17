library(DESeq2)
install.packages("tidyverse")
BiocManager::install("apeglm")
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
res <- results(dds)
summary(res)

res0.01 <- results(dds, alpha = 0.05)
summary(res0.01)
res

plotMA(res)

vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = "condition")

resultsNames(dds)
resLFC <- lfcShrink(dds, coef = "condition_snf2_vs_WT", type = "apeglm")

plotMA(resLFC)

plot(
  res$log2FoldChange,
  resLFC$log2FoldChange,
  xlab = "Unshrunken LFC",
  ylab = "Shrunken LFC"
)
abline(0, 1, col = "red")


resdf <- as.data.frame(resLFC)
resdf$sig <- !is.na(resdf$padj) & 
             resdf$padj < 0.05 & 
             abs(resdf$log2FoldChange) > 1

ggplot(resdf, aes(log2FoldChange, -log10(pvalue), colour = sig)) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_colour_manual(values = c("grey70", "#c0392b")) +
  theme_minimal() + labs(colour = "padj<0.05 & |LFC|>1")

write.csv(as.data.frame(resLFC[order(resLFC$padj), ]),
          "../results/de_snf2_vs_WT.csv")

resLFC["YOR290C", ]
