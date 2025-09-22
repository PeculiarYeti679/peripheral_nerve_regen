# 4 principal component analysis

#load libraries
library(ggplot2)
library(dplyr)
library(pheatmap)
library(tibble)

#load matrix from the file saved in 1_preprocessing.R
expression_matrix <- readRDS("expression_matrix.rds")

# clean up sample names before PCA and heatmap
clean_names <- sub("^GSM\\d+_", "", colnames(expression_matrix))  # remove GSM ID
clean_names <- sub("\\.txt$", "", clean_names)                    # remove .txt
colnames(expression_matrix) <- clean_names

# run pca from a transposed matrix
pca <- prcomp(t(expression_matrix), scale. = TRUE)

#convert to a data frame and then merge it with the phenodata created earlier
pca_df <- as.data.frame(pca$x) 
pca_df$Sample <- rownames(pca_df)

#load in the phenodata
phenodata <- read.csv("phenodata.csv", stringsAsFactors = FALSE)

#also clean Sample column in phenodata to match
phenodata$Sample <- sub("^GSM\\d+_", "", phenodata$Sample)
phenodata$Sample <- sub("\\.txt$", "", phenodata$Sample)

#merge the data frame from above with the loaded phenodata
pca_df <- left_join(pca_df, phenodata, by = "Sample")
head(pca_df)
#plot using the the tissue for color
ggplot(pca_df, aes(x = PC1, y = PC2, color = Tissue, label = Timepoint)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(size = 3, vjust = 1.6, color = "black") +
  labs(title = "PCA of Samples Colored by Tissue", x = "PC1", y = "PC2") +
  theme_minimal()

#plot using the timepoints to color
ggplot(pca_df, aes(x = PC1, y = PC2, color = Timepoint, label = Tissue)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(size = 3, vjust = 1.6, color = "black") +
  labs(title = "PCA of Samples Colored by Timepoint", x = "PC1", y = "PC2") +
  theme_minimal()

#show correlation between samples using a matrix
sample_cor <- cor(expression_matrix, method = "pearson")
print(sample_cor)
#annotate the matrix using phenodata
annotation <- phenodata %>% column_to_rownames("Sample")
head(annotation)

#create the heatmap using the correlation matrix and the annotated columns
pheatmap(sample_cor,
         annotation_col = annotation[, c("Tissue", "Timepoint")],
         main = "Sample-to-Sample Correlation",
         clustering_method = "ward.D2",
         border_color = NA)

#export data from above to use later and in the frontend
write.csv(expression_matrix, "expression_matrix.csv")
write.csv(phenodata, "phenodata.csv", row.names = FALSE)
saveRDS(expression_matrix, "expression_matrix.rds")
saveRDS(phenodata, "phenodata.rds")
ggsave("PCA_plot_by_tissue.png", width = 8, height = 6, dpi = 300)
pheatmap(sample_cor,
         annotation_col = annotation[, c("Tissue", "Timepoint")],
         filename = "sample_correlation_heatmap.png",
         width = 8,
         height = 8)



#export data for use in frontend visualizations
write.csv(expression_matrix, "frontend_expression_matrix.csv")
write.csv(pca_df, "frontend_pca_coordinates.csv", row.names = FALSE)
write.csv(sample_cor, "frontend_sample_correlation.csv")
write.csv(phenodata, "frontend_phenodata.csv", row.names = FALSE)

#l: also export as JSON
library(jsonlite)
write_json(as.data.frame(expression_matrix), "frontend_expression_matrix.json", pretty = TRUE)
write_json(pca_df, "frontend_pca_coordinates.json", pretty = TRUE)
write_json(sample_cor, "frontend_sample_correlation.json", pretty = TRUE)
write_json(phenodata, "frontend_phenodata.json", pretty = TRUE)

