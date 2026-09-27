#!/usr/bin/env Rscript

library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(org.Sc.sgd.db)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep("^--file=", args)
  if (length(match) > 0) return(dirname(normalizePath(sub("^--file=", "", args[match]))))
  return(getwd())
}
theme_file <- file.path(get_script_dir(), "plot_theme.R")
if (file.exists(theme_file)) source(theme_file)
if (exists("custom_theme")) theme_set(custom_theme)

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
  
  ego_simple <- tryCatch(
    simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min),
    error = function(e) ego
  )

  p1 <- dotplot(ego_simple, showCategory = 10, label_format = 40) + 
    ggtitle("ORA: GO Biological Process") +
    theme(axis.text.y = element_text(size = 10, face = "bold"))
  ggsave("enrichment_output/ORA_dotplot.png", plot = p1, width = 9, height = 6.5)
  
  p2 <- cnetplot(ego_simple, showCategory = 5)
  ggsave("enrichment_output/ORA_cnetplot.png", plot = p2, width = 9, height = 7)
  
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
    quit(save = "no", status = 1)
  }, finally = {
    cat("Enrichment analysis has been executed.\n")
  })
  