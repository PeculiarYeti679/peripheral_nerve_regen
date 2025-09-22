# 3 phenodata generation

# load matrix from the file saved in 1_preprocessing.R
expression_matrix <- readRDS("expression_matrix.rds")

# get sample names from the expression matrix
sample_names <- colnames(expression_matrix)
# check the information
head(sample_names)

# parse the data using regex
phenodata <- data.frame(
  Sample    = sample_names,  # Always matches colnames
  GSM       = paste0("Sample_", seq_along(sample_names)),  # Placeholder GSM
  Timepoint = sub("_.*", "", sample_names),                 # Extract timepoint
  Tissue    = sub(".*_(DRG|SN)_.*", "\\1", sample_names),   # Extract DRG or SN
  Replicate = sub(".*_(\\d+)$", "\\1", sample_names),       # Extract replicate
  stringsAsFactors = FALSE
)

# Convert Timepoint into a factor with chronological order
phenodata$Timepoint <- factor(
  phenodata$Timepoint,
  levels = c("0h", "1d", "4d", "7d", "14d"),
  ordered = TRUE
)

# Convert Tissue to a factor (DRG, SN)
phenodata$Tissue <- factor(phenodata$Tissue, levels = c("DRG", "SN"))

# validate alignment
stopifnot(identical(phenodata$Sample, colnames(expression_matrix)))

# check the phenodata
head(phenodata)

# write the data to a CSV
write.csv(phenodata, "phenodata.csv", row.names = FALSE)

# check balance
table(phenodata$Timepoint, phenodata$Tissue)

