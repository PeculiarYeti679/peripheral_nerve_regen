# 6 Differentially Expressed Genes (DEGs)  
# -----------------------------
# Load libraries
library(tidyverse)
library(affy)        
library(limma)
library(clusterProfiler)
library(org.Rn.eg.db)
library(pheatmap)
library(ggplot2)


#current dataset uses Agilent 1-color arrays

raw_dir <- "./GSE30165_raw/CEL/"
files <-  list.files(raw_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)