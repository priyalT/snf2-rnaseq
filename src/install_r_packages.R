# Installs the R packages used by count_matrix.R and run_deseq2.R.
#
#   Rscript src/install_r_packages.R
#
# Versions that produced the results in this repo (see README):
#   R 4.6.0 (2026-04-24)   DESeq2 1.52.0   apeglm 1.34.0

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

BiocManager::install(c("DESeq2", "apeglm"), ask = FALSE, update = FALSE)

install.packages(
  c("tidyverse", "readr", "dplyr", "glue", "ggplot2"),
  repos = "https://cloud.r-project.org"
)

pkgs <- c("DESeq2", "apeglm", "tidyverse", "readr", "dplyr", "glue", "ggplot2")
versions <- vapply(pkgs, function(p) as.character(utils::packageVersion(p)), character(1))
writeLines(
  c(paste("R", getRversion()),
    paste("Bioconductor", BiocManager::version()),
    paste(pkgs, versions)),
  "../results/r_package_versions.txt"
)
print(data.frame(package = pkgs, version = versions))
