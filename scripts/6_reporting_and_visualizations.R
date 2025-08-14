#7 reporting and visualizations
# ===========================================
# Phase 5: Final Analysis, Visualization & Export
# ===========================================

# -----------------------------
# 1. Load Libraries
# -----------------------------
library(tidyverse)
library(limma)
library(pheatmap)
library(clusterProfiler)
library(org.Rn.eg.db)
library(ggplot2)

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