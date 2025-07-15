#quality check (qc) visualization

#load libraries
library(readr)
library(ggplot2)
library(dplyr)
library(pheatmap)

#load matrix from the file saved in 1_preprocessing.R
expression_matrix <- readRDS("expression_matrix.rds")

#checking data(qc)
dim(expression_matrix)     
head(expression_matrix[, 1:5]) 
boxplot(expression_matrix[, 1:5]) 
#summary of entire data set 
summary(as.vector(as.matrix(expression_matrix)))
#count all NA values per sample
colSums(is.na(expression_matrix))

#count all NA values in a probe
rowSums(is.na(expression_matrix)) |> summary()

#boxplot of expression values for each sample
boxplot(expression_matrix,
        las = 2,
        main = "Expression Distribution by Sample",
        ylab = "Log2 Expression (gProcessedSignal)",
        col = "lightblue",
        outline = FALSE)

#PCA - Principal Component Analysis
#Looking for outliers from the main cluster
pca <- prcomp(t(expression_matrix), scale. = TRUE)

# Prepare a data frame for plotting
pca_df <- as.data.frame(pca$x)
pca_df$Sample <- rownames(pca_df)

ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(color = "darkred") +
  geom_text(size = 3, vjust = 1.5) +
  labs(title = "PCA of Gene Expression", x = "PC1", y = "PC2")

