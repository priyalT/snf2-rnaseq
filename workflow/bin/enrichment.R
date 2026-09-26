#!/usr/bin/env Rscript

library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(org.Sc.sgd.db)


args <- commandArgs(trailingOnly = TRUE)

tryCatch({
  
  deseq <- read.csv(args[1], header=TRUE, row.names="X")

  
  #ORA
  sig_genes <- rownames(subset(deseq, padj < 0.05 & abs(log2FoldChange) > 1))
  universe <- rownames(deseq)
  ego <- enrichGO(gene          = sig_genes,
                  universe      = universe,
                  OrgDb         = org.Sc.sgd.db,
                  keyType       = "ORF",
                  ont           = "BP",          
                  pAdjustMethod = "BH",
                  qvalueCutoff  = 0.05,
                  readable      = FALSE)          
  
  dir.create("enrichment_output", showWarnings = FALSE)
  
  p1 <- dotplot(ego, showCategory = 12) + ggtitle("ORA — GO biological process")
  ggsave("enrichment_output/ORA_dotplot.png", plot = p1, width = 8, height = 6)
  
  p2 <- cnetplot(ego, showCategory = 5)
  ggsave("enrichment_output/ORA_cnetplot.png", plot = p2, width = 8, height = 6)
  
  write.csv(as.data.frame(ego), "enrichment_output/ORA_results.csv")
  
  #GSEA
  ranks <- deseq$log2FoldChange
  names(ranks) <- rownames(deseq)
  ranks <- sort(ranks[!is.na(ranks)], decreasing = TRUE)
  set.seed(42)                               
  gse <- gseGO(geneList      = ranks,
               OrgDb         = org.Sc.sgd.db,
               keyType       = "ORF",
               ont           = "BP",
               pAdjustMethod = "BH",
               verbose       = FALSE)
  gdf <- as.data.frame(gse)
  gdf <- gdf[order(-abs(gdf$NES)), ]
  
  
  top_id <- gdf$ID[which.max(gdf$NES)]
  p3 <- gseaplot2(gse, geneSetID = top_id, title = gdf$Description[gdf$ID == top_id])
  ggsave("enrichment_output/GSEA_top_pathway.png", plot = p3, width = 8, height = 6)
  
  write.csv(gdf, "enrichment_output/GSEA_results.csv")
},
  error = function(e) {
    cat("An error occured: ", conditionMessage(e), "\n")
    NA
  }, finally = {
    cat("Enrichment analysis has been executed.\n")
  })
  