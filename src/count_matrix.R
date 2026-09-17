library(readr)
library(dplyr)
library(glue)
library(DESeq2)
library(ggplot2)

coldata <- read.csv("../data/samplesheet.csv")
samples <- coldata$sample

for (sample in samples) {
  counts  <- as.data.frame(
    read_delim(glue("../results/alignment/star/{sample}_ReadsPerGene.out.tab"),
               delim = "\t", col_names = FALSE))
  counts_df <- as.data.frame(counts[, -1])
  rownames(counts_df) <- counts$X1
  counts_df_final <- counts_df %>%
    select(1)
  colnames(counts_df_final)[colnames(counts_df_final) == 'X2'] <- glue("{sample}")
  assign(paste0("counts_", sample), counts_df_final)
}

merged_df <- Reduce(
  function(x, y) {
    
    z <- merge(x, y, by = "row.names", all = TRUE)
    
    rownames(z) <- z$Row.names
    z$Row.names <- NULL
    
    return(z)
  },
  list(
    counts_ERR458495,
    counts_ERR458502,
    counts_ERR458509,
    counts_ERR458517,
    counts_ERR458528,
    counts_ERR458535,
    counts_ERR458552,
    counts_ERR458882,
    counts_ERR458887,
    counts_ERR458905,
    counts_ERR458906,
    counts_ERR458921
  )
)

class(merged_df)
write.csv(merged_df, 
          "../results/countmatrix.csv")