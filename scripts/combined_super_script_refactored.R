# ==============================
# Load Required Libraries
# ==============================
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(jsonlite)
library(limma)
library(org.Rn.eg.db)
library(pheatmap)
library(purrr)
library(readr)
library(tibble)
library(tidyverse)
library(tools)



# =============================================
# START OF 1_preprocessing.R
# =============================================

#preprocess data

#download libraries
library(readr) #read in tab delimited txt files
library(purrr) #functional programming; reduce()
library(dplyr) #data manipulation; summarise 
library(tools) #file name manipulation 

#get all files in the untarred CEL folder with .txt.gz file type
files <- list.files(".", pattern = "*.txt.gz", full.names = TRUE)

#check for that amount of files in the folder (30) in this case
length(files)

#read the files in and move to dataframes
data_list <- lapply(files, function(f) {
  #this is the sample name
  sample_name <- file_path_sans_ext(basename(f))
  df <- read_tsv(
    f,
    #skip the first 9 lines that do no contain actual data 
    skip = 9,
    show_col_types = FALSE,
    col_types = cols_only(
      ProbeName = col_character(),
      gProcessedSignal = col_double()
    )
  )
  
  # collapse duplicate ProbeIDs by averaging their expression
  df_clean <- df %>%
    group_by(ProbeName) %>%
    summarise(!!sample_name := mean(gProcessedSignal, na.rm = TRUE), .groups = "drop")
  
  #rename the "ProbeName" to the "ProbeID" so that all samples can be consistent
  colnames(df_clean)[1] <- "ProbeID"
  return(df_clean)
})
print(data_list)

# merge all samples into one matrix
expression_matrix <- reduce(data_list, ~ inner_join(.x, .y, by = "ProbeID"))



# move ProbeID to rownames
rownames(expression_matrix) <- expression_matrix$ProbeID
expression_matrix <- expression_matrix[, -1]

#save the files to reuse later
saveRDS(expression_matrix, file = "expression_matrix.rds")





# =============================================
# START OF 2_qc_visualization.R
# =============================================

#quality check (qc) visualization

#load libraries
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





# =============================================
# START OF 3_phenodata_generation.R
# =============================================

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



# =============================================
# START OF 4_principal_component_analysis.R
# =============================================

# 4 principal component analysis

#load libraries
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
write_json(as.data.frame(expression_matrix), "frontend_expression_matrix.json", pretty = TRUE)
write_json(pca_df, "frontend_pca_coordinates.json", pretty = TRUE)
write_json(sample_cor, "frontend_sample_correlation.json", pretty = TRUE)
write_json(phenodata, "frontend_phenodata.json", pretty = TRUE



# =============================================
# START OF 5_differential_expression.R
# =============================================

# 5 differential expression 
#load in generated data
expression_matrix <- readRDS("expression_matrix.rds")
phenodata <- read.csv("phenodata.csv", stringsAsFactors = FALSE)

#normalize to log2 for limma
expression_matrix <- log2(expression_matrix + 1)
hist(as.numeric(as.matrix(expression_matrix)),
     main = "Histogram of Log2 Expression Values",
     xlab = "Log2 Expression",
     col = "skyblue",
     breaks = 50)


#clean the names of the phenodata 
phenodata$Sample <- sub("^GSM\\d+_", "", phenodata$Sample)
phenodata$Sample <- sub("\\.txt$", "", phenodata$Sample)

#reorder the phenodata to match the expression matrix columns
phenodata <- phenodata[match(colnames(expression_matrix), phenodata$Sample), ]

#create factor - for categorical variable tissue in the phenodata (DRG - SN)
group <- factor(phenodata$Tissue)

#create matrix
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

#linear model fit
fit <- lmFit(expression_matrix, design)

#looking at DRG vs SN
contrast.matrix <- makeContrasts(DRG_vs_SN = DRG - SN, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

topTable_DRG_vs_SN <- topTable(fit2, coef = "DRG_vs_SN", number = Inf, adjust.method = "fdr")
head(topTable_DRG_vs_SN)

#save results
write.csv(topTable_DRG_vs_SN, "limma_DRG_vs_SN_results.csv")

 #filter genes -  FDR < 0.05 and absolute logFC > 1
de_genes <- topTable_DRG_vs_SN %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1)
#save filtered results
write.csv(de_genes, "filtered_DE_genes.csv", row.names = TRUE)

#get top 10 by lowest adjusted p-value
top10_genes <- topTable_DRG_vs_SN %>%
  arrange(adj.P.Val) %>%
  head(10)
#save top 10 genes
write.csv(top10_genes, "top10_DE_genes.csv", row.names = TRUE)


# Add significance categories for coloring
topTable_DRG_vs_SN$Significant <- with(topTable_DRG_vs_SN,
                                       ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "Significant", "Not Significant")
)

# Add a column to flag top 10 genes
topTable_DRG_vs_SN$Top10 <- rownames(topTable_DRG_vs_SN) %in% rownames(top10_genes)

# Plot
ggplot(topTable_DRG_vs_SN, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Significant), alpha = 0.6) +
  geom_point(data = subset(topTable_DRG_vs_SN, Top10 == TRUE),
             color = "red", size = 2.5) +
  geom_text(data = subset(topTable_DRG_vs_SN, Top10 == TRUE),
            aes(label = rownames(subset(topTable_DRG_vs_SN, Top10 == TRUE))),
            vjust = -1, size = 3.5) +
  scale_color_manual(values = c("grey", "blue")) +
  labs(title = "Volcano Plot: DRG vs SN",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-value") +
  theme_minimal()

# Save volcano plot
ggsave("volcano_plot_DRG_vs_SN.png", width = 8, height = 6, dpi = 300)

# Export full table
write_json(topTable_DRG_vs_SN, "DE_results_all.json", pretty = TRUE)

# Optional: Export filtered significant genes
write_json(de_genes, "DE_results_filtered.json", pretty = TRUE)

# Optional: Export top 10 separately
write_json(top10_genes, "DE_results_top10.json", pretty = TRUE)

topTable_DRG_vs_SN$Gene <- rownames(topTable_DRG_vs_SN)
topTable_json <- topTable_DRG_vs_SN %>% select(Gene, everything())
write_json(topTable_json, "DE_results_all_labeled.json", pretty = TRUE)



# =============================================
# START OF 6_differential_expression_timepoint.R
# =============================================

# ===========================================
# 6. Differentially Expressed Genes (DEGs)
# Phase 4: Temporal and Functional Pattern Analysis
# ===========================================

# -----------------------------
# Load libraries
# -----------------------------
# 1. Prepare Expression Data
# -----------------------------

# Remove ".txt" extension from column names
colnames(expression_matrix) <- sub("\\.txt$", "", colnames(expression_matrix))

# Create phenodata from column names
phenodata <- data.frame(Sample = colnames(expression_matrix)) %>%
  separate(Sample, into = c("GSM", "Timepoint", "Tissue", "Replicate"), sep = "_")

# -----------------------------
# 2. Design Matrix for Limma
# -----------------------------
phenodata$Tissue <- factor(phenodata$Tissue, levels = c("DRG", "SN"))
design <- model.matrix(~ 0 + Tissue, data = phenodata)
colnames(design) <- c("DRG", "SN")
print(design)

# -----------------------------
# 3. Linear Modeling with Limma
# -----------------------------
fit <- lmFit(expression_matrix, design)
contrast.matrix <- makeContrasts(DRG_vs_SN = DRG - SN, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

topTable_DRG_vs_SN <- topTable(fit2, coef = "DRG_vs_SN", number = Inf)
head(topTable_DRG_vs_SN)
cat("Number of significant DEGs (adj.P.Val < 0.05):", 
    sum(topTable_DRG_vs_SN$adj.P.Val < 0.05), "\n")

# -----------------------------
# 4. Filter DEGs
# -----------------------------
deg_filtered <- topTable_DRG_vs_SN %>%
  filter(adj.P.Val < 0.05) %>%
  arrange(adj.P.Val)

deg_probes <- deg_filtered$ID
cat("Number of DEG probes:", length(deg_probes), "\n")

# Subset expression matrix for DEGs
expression_deg <- expression_matrix[deg_probes, ]
dim(expression_deg)

# -----------------------------
# 5. Collapse Duplicates
# -----------------------------
expression_deg_unique <- expression_deg %>%
  as.data.frame() %>%
  rownames_to_column("ProbeID") %>%
  group_by(ProbeID) %>%
  summarize(across(everything(), mean), .groups = "drop") %>%
  column_to_rownames("ProbeID")

cat("Duplicated probes after collapsing:", 
    sum(duplicated(rownames(expression_deg_unique))), "\n")
dim(expression_deg_unique)

# -----------------------------
# 6. Prepare Long Data for Clustering
# -----------------------------
expr_transposed <- as.data.frame(t(expression_deg_unique))
expr_transposed$Sample <- rownames(expr_transposed)

# Ensure phenodata includes "Sample"
phenodata <- phenodata %>%
  mutate(Sample = paste(GSM, Timepoint, Tissue, Replicate, sep = "_"))

# Join with metadata
expr_long <- expr_transposed %>%
  left_join(phenodata, by = "Sample")

# Pivot longer (Gene/Expression columns)
gene_cols <- setdiff(colnames(expr_long), c("Sample", "GSM", "Timepoint", "Tissue", "Replicate"))
expr_long <- expr_long %>%
  pivot_longer(
    cols = all_of(gene_cols),
    names_to = "Gene",
    values_to = "Expression"
  )

# -----------------------------
# 7. Preview Final Data
# -----------------------------
head(expr_long)

dim(expr_long)
head(expr_long)

# Average expression for each Gene-Timepoint-Tissue
expr_timepoint <- expr_long %>%
  group_by(Gene, Timepoint, Tissue) %>%
  summarise(Expression = mean(Expression), .groups = "drop")

head(expr_timepoint)

expr_wide <- expr_timepoint %>%
  unite("Time_Tissue", Timepoint, Tissue, sep = "_") %>%
  pivot_wider(names_from = Time_Tissue, values_from = Expression)

# Check structure
dim(expr_wide)
head(expr_wide)

# Prepare matrix (remove Gene column)
gene_matrix <- as.matrix(expr_wide[,-1])
rownames(gene_matrix) <- expr_wide$Gene

# Hierarchical clustering heatmap
pheatmap(
  gene_matrix,
  scale = "row",            # Z-score per gene
  clustering_method = "complete",
  show_rownames = FALSE,
  show_colnames = TRUE
)

head(raw_data$genes)
# Create a mapping of probe IDs to gene names
probe_to_gene <- raw_data$genes %>%
  select(ProbeName, GeneName) %>%
  distinct()

# Merge gene names into expression matrix
expression_deg_annotated <- expression_deg_unique %>%
  rownames_to_column("ProbeID") %>%
  left_join(probe_to_gene, by = c("ProbeID" = "ProbeName"))

head(expression_deg_annotated)

# Convert probe IDs to Entrez IDs
deg_entrez <- mapIds(org.Rn.eg.db,
                     keys = expression_deg_annotated$GeneName,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")

# Add Entrez IDs to the dataframe
expression_deg_annotated$EntrezID <- deg_entrez
head(expression_deg_annotated)

deg_entrez <- na.omit(deg_entrez)

# GO enrichment
ego <- enrichGO(gene         = deg_entrez,
                OrgDb        = org.Rn.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",
                pAdjustMethod= "BH",
                qvalueCutoff = 0.05,
                readable     = TRUE)

# KEGG enrichment
ekegg <- enrichKEGG(gene     = deg_entrez,
                    organism = "rno",
                    pAdjustMethod = "BH",
                    qvalueCutoff  = 0.05)

pheatmap(expression_deg_unique, scale = "row")
dotplot(ego, showCategory = 20)
dotplot(ekegg, showCategory = 20)
write.csv(as.data.frame(ego), "GO_enrichment_results.csv")
write.csv(as.data.frame(ekegg), "KEGG_enrichment_results.csv")

# Final DEG table
deg_final <- topTable_DRG_vs_SN %>%
  rownames_to_column("ProbeID") %>%
  left_join(probe_to_gene, by = c("ProbeID" = "ProbeName")) %>%
  mutate(EntrezID = mapIds(org.Rn.eg.db,
                           keys = GeneName,
                           column = "ENTREZID",
                           keytype = "SYMBOL",
                           multiVals = "first"))

write.csv(deg_final, "DEG_final_annotated.csv", row.names = FALSE)
# ------------------------------------------
# 8. Final Visualization & Export
# ------------------------------------------

# Volcano plot
volcano_data <- topTable_DRG_vs_SN %>%
  rownames_to_column("ProbeID") %>%
  mutate(Significant = adj.P.Val < 0.05)

ggplot(volcano_data, aes(x = logFC, y = -log10(P.Value), color = Significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: DRG vs SN", x = "log2 Fold Change", y = "-log10(P.Value)")
ggsave("volcano_plot_DRG_vs_SN.png", width = 6, height = 5)

# PCA
pca <- prcomp(t(expression_matrix), scale. = TRUE)
pca_df <- data.frame(pca$x, phenodata)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Tissue, shape = Timepoint)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA of Samples", x = "PC1", y = "PC2")
ggsave("PCA_plot.png", width = 6, height = 5)

# GO & KEGG barplots
barplot(ego, showCategory = 20, title = "Top 20 GO Biological Processes")
ggsave("GO_enrichment_barplot.png", width = 7, height = 5)

barplot(ekegg, showCategory = 20, title = "Top 20 KEGG Pathways")
ggsave("KEGG_enrichment_barplot.png", width = 7, height = 5)

# Top 50 DEG heatmap
top50_probes <- deg_filtered$ID[1:50]
top50_matrix <- expression_matrix[top50_probes, ]

pheatmap(
  top50_matrix,
  scale = "row",
  annotation_col = phenodata[, c("Timepoint", "Tissue")],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE
)
ggsave("heatmap_top50_DEGs.png", width = 8, height = 6)

# Export CSVs
write.csv(expr_long, "Expression_Long_Format.csv", row.names = FALSE)
write.csv(expr_timepoint, "Expression_Timepoint_Averages.csv", row.names = FALSE)
write.csv(deg_final, "DEG_final_annotated.csv", row.names = FALSE)
write.csv(as.data.frame(ego), "GO_enrichment_results.csv")
write.csv(as.data.frame(ekegg), "KEGG_enrichment_results.csv")



# =============================================
# START OF 6_reporting_and_visualizations.R
# =============================================

#7 reporting and visualizations
# ===========================================
# Phase 5: Final Analysis, Visualization & Export
# ===========================================

# -----------------------------
# 1. Load Libraries
# -----------------------------
# -----------------------------
# 2. Load Processed Data
# -----------------------------
# Replace these paths with your saved files
expression_matrix <- readRDS("agilent_expression_matrix.rds")
topTable_DRG_vs_SN <- readRDS("topTable_DRG_vs_SN.rds")
phenodata <- read.csv("./csv/phenodata.csv", stringsAsFactors = FALSE)
ego <- readRDS("GO_enrichment.rds")
ekegg <- readRDS("KEGG_enrichment.rds")

# -----------------------------
# 3. DEG Table with Annotations
# -----------------------------
probe_to_gene <- readRDS("probe_to_gene.rds")  # Ensure you saved this mapping

deg_final <- topTable_DRG_vs_SN %>%
  rownames_to_column("ProbeID") %>%
  left_join(probe_to_gene, by = c("ProbeID" = "ProbeName")) %>%
  mutate(EntrezID = mapIds(org.Rn.eg.db,
                           keys = GeneName,
                           column = "ENTREZID",
                           keytype = "SYMBOL",
                           multiVals = "first"))

write.csv(deg_final, "results/DEG_final_annotated.csv", row.names = FALSE)

# -----------------------------
# 4. Volcano Plot
# -----------------------------
volcano_data <- topTable_DRG_vs_SN %>%
  rownames_to_column("ProbeID") %>%
  mutate(Significant = adj.P.Val < 0.05)

ggplot(volcano_data, aes(x = logFC, y = -log10(P.Value), color = Significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: DRG vs SN", x = "log2 Fold Change", y = "-log10(P.Value)")

ggsave("results/volcano_plot_DRG_vs_SN.png", width = 6, height = 5)

# -----------------------------
# 5. PCA of Samples
# -----------------------------
pca <- prcomp(t(expression_matrix), scale. = TRUE)
pca_df <- data.frame(pca$x, phenodata)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Tissue, shape = Timepoint)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA of Samples", x = "PC1", y = "PC2")

ggsave("results/PCA_plot.png", width = 6, height = 5)

# -----------------------------
# 6. GO and KEGG Enrichment
# -----------------------------
barplot(ego, showCategory = 20, title = "Top 20 GO Biological Processes")
ggsave("results/GO_enrichment_barplot.png", width = 7, height = 5)

barplot(ekegg, showCategory = 20, title = "Top 20 KEGG Pathways")
ggsave("results/KEGG_enrichment_barplot.png", width = 7, height = 5)

write.csv(as.data.frame(ego), "results/GO_enrichment_results.csv")
write.csv(as.data.frame(ekegg), "results/KEGG_enrichment_results.csv")

# -----------------------------
# 7. Heatmap of Top 50 DEGs
# -----------------------------
deg_filtered <- topTable_DRG_vs_SN %>%
  filter(adj.P.Val < 0.05) %>%
  arrange(adj.P.Val)
top50_probes <- deg_filtered$ID[1:50]

top50_matrix <- expression_matrix[top50_probes, ]

pheatmap(
  top50_matrix,
  scale = "row",
  annotation_col = phenodata[, c("Timepoint", "Tissue")],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE
)

ggsave("results/heatmap_top50_DEGs.png", width = 8, height = 6)

# -----------------------------
# 8. Save Additional Outputs
# -----------------------------
write.csv(expression_matrix, "results/Normalized_Expression_Matrix.csv")

