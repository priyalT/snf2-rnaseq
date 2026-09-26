#!/usr/bin/env Rscript

library(DESeq2)
library(dplyr)
library(readr)
library(here)
library(ggplot2)
library(apeglm)
args <- commandArgs(trailingOnly = TRUE)
counts_data <- read.csv("../results/counts/countmatrix.csv", header=TRUE, row.names="X")
