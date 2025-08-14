
#create directories to store data
dirs <- c(
  "data/raw",        # downloads from GEO
  "data/processed",  # intermediate .rds objects
  "data/meta",       # phenodata/annotation tables
  "outputs/qc",      # QC plots
  "outputs/tables",  # CSV tables for DEGs, etc.
  "outputs/figures"  # final figures
)
invisible(lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE))

#load libraries for fetching data, reading the microarray, normalziing
#and data manipulation
library(GEOquery)
library(limma)
library(tidyverse)

#get data
gse_id <- "GSE30165"

gse_list <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE)

length(gse_list) 
eset <- gse_list[[1]]
#phenotype data (pd)
pd <- Biobase::pData(eset)
write.csv(pd, file = "data/meta/phenodata_raw.csv", row.names = TRUE)

#check pd 
dim(pd)
head(colnames(pd), 20)
head(pd[, 1:6])

#create phenodata table 
library(tidyverse)
library(tibble)
pdTable <- read.csv(
  "data/meta/phenodata_raw.csv",
  check.names = FALSE,
  row.names = 1       
)
colnames(pdTable)
keep <- pd %>%
  select(geo_accession, title, source_name_ch1, starts_with("characteristics_ch1"))
head(keep)

#make long strings easier to read via parsing
parse_tissue <- function(x) {
  x <- as.character(x)
  out <- str_match(x, "^\\s*([A-Za-z]+)\\s*_")[,2]
  str_to_upper(out)
}

parse_timepoint <- function(x) {
  x <- as.character(x)
  m <- str_match(x, "([0-9]+(?:\\.[0-9]+)?)\\s*([hd])\\b")  # captures number + unit
  ifelse(is.na(m[,1]), NA, paste0(m[,2], m[,3]))           # "0d", "0.5h", "3h"
}

time_to_hours <- function(tstr) {
  m <- str_match(tstr, "^\\s*([0-9]+(?:\\.[0-9]+)?)([hd])\\s*$")
  val <- suppressWarnings(as.numeric(m[,2]))
  unit <- m[,3]
  ifelse(is.na(val), NA_real_,
         ifelse(unit == "d", val * 24, val))
}

parse_injury_phrase <- function(x) {
  x <- as.character(x)
  # everything between first and last underscore
  out <- str_match(x, "^[^_]+_(.*)_[^_]+$")[,2]
  ifelse(is.na(out), NA, str_squish(out))
}

# Build phenodata ----
phenodata <- keep %>%
  mutate(
    Sample       = geo_accession,
    raw_source   = source_name_ch1,
    tissue       = parse_tissue(raw_source),
    time         = parse_timepoint(raw_source),
    time_hours   = time_to_hours(time),
    injury_label = parse_injury_phrase(raw_source)
  ) %>%
  transmute(Sample, tissue, time, time_hours, injury_label, raw_source)

#check data before saving to CSV
print(phenodata, 8)
table(tissue = phenodata$tissue, useNA = "ifany")
table(time   = phenodata$time,   useNA = "ifany")

#save the processed data
write.csv(phenodata, "data/meta/phenodata.csv", row.names = FALSE)

#get expression set
expr_processed <- Biobase::exprs(eset)
dim(expr_processed)
expr_processed[1:3, 1:3]

pheno <- read.csv("data/meta/phenodata.csv")
expr_processed <- expr_processed[, pheno$Sample, drop = FALSE]
print(expr_processed)

saveRDS(expr_processed, file = "data/processed/expression_matrix_processed.rds")


#qc the data from above

#load in plotting library
library(ggplot2)
pheno <- read.csv("data/meta/phenodata.csv")

exprQC <- readRDS("data/processed/expression_matrix_processed.rds")

boxplot(exprQC, outline = FALSE, main = "Expression Distributions", las = 2)

expr_t <- t(exprQC) 
pca <- prcomp(expr_t, scale. = TRUE)
pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
) %>%
  left_join(pheno, by = "Sample")
ggplot(pca_df, aes(PC1, PC2, label = Sample, shape = tissue)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(min.segment.length = 0) 
library(ggrepel)
pca_df <- pca_df %>%
  mutate(label_flag = ifelse(abs(PC1) > 20 | abs(PC2) > 20, Sample, NA))

ggplot(pca_df, aes(PC1, PC2, shape = tissue)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = label_flag), min.segment.length = 0) +
  theme_minimal()
ggplot(pca_df, aes(PC1, PC2, color = tissue, shape = tissue)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Expression Data")
