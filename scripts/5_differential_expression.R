# 5 differential expression 
library(limma)
library(dplyr)

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


library(ggplot2)

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

library(jsonlite)

# Export full table
write_json(topTable_DRG_vs_SN, "DE_results_all.json", pretty = TRUE)

# Optional: Export filtered significant genes
write_json(de_genes, "DE_results_filtered.json", pretty = TRUE)

# Optional: Export top 10 separately
write_json(top10_genes, "DE_results_top10.json", pretty = TRUE)

topTable_DRG_vs_SN$Gene <- rownames(topTable_DRG_vs_SN)
topTable_json <- topTable_DRG_vs_SN %>% select(Gene, everything())
write_json(topTable_json, "DE_results_all_labeled.json", pretty = TRUE)