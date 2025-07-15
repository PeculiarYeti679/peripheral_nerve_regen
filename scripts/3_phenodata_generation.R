#3 phenodata generation

#load matrix from the file saved in 1_preprocessing.R
expression_matrix <- readRDS("expression_matrix.rds")

#get sample names from the expression matrix
sample_names <- colnames(expression_matrix)
#check the information 
head(sample_names)

#parse the data using regex
#drg = dorsal root ganglia 
#sn = sciatic nerve
phenodata <- data.frame(
  Sample = sample_names,
  GSM = sub("_.*", "", sample_names),  # Everything before first underscore
  Timepoint = sub(".*_(\\d+h|\\d+d)_.*", "\\1", sample_names),  # 0h, 1h, 14d etc.
  Tissue = sub(".*_(DRG|SN)_.*", "\\1", sample_names),  # DRG or SN
  Replicate = sub(".*_(\\d+)\\.txt", "\\1", sample_names),  # Final digit before .txt
  stringsAsFactors = FALSE
)
#check the phenodata
head(phenodata)

#write the data to a csv
write.csv(phenodata, "phenodata.csv", row.names = FALSE)

#check balance
table(phenodata$Timepoint, phenodata$Tissue)