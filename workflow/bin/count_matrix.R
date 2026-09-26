#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(glue)
library(DESeq2)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
coldata <- read.csv(args[1])
samples <- coldata$sample

for (sample in samples) {
  counts  <- as.data.frame(
    read_delim(glue("{sample}_trimmed_ReadsPerGene.out.tab"),
               delim = "\t", col_names = FALSE))
  counts_df <- as.data.frame(counts[, -1])
  rownames(counts_df) <- counts$X1
  counts_df_final <- counts_df %>%
    select(1)
  colnames(counts_df_final)[colnames(counts_df_final) == 'X2'] <- glue("{sample}")
  assign(paste0("counts_", sample), counts_df_final)
}
counts_list <- lapply(samples, function(s) get(paste0("counts_", s)))

merged_df <- Reduce(
  function(x, y) {
    
    z <- merge(x, y, by = "row.names", all = TRUE)
    
    rownames(z) <- z$Row.names
    z$Row.names <- NULL
    
    return(z)
  },
  counts_list
)

merged_df <- merged_df[!row.names(merged_df) %in% c("N_ambiguous", "N_multimapping",
                                    "N_noFeature", "N_unmapped"), , drop = FALSE]


write.csv(merged_df, "countmatrix.csv")