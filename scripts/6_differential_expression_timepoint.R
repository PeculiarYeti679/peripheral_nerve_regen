# ===========================================
# 6. Differentially Expressed Genes (DEGs)
# ===========================================

# -----------------------------
# Load required libraries
# -----------------------------
library(limma)
library(dplyr)
library(ggplot2)

# -----------------------------
# 1. Load Data
# -----------------------------
# Load expression matrix (from 1_preprocessing.R)
expression_matrix <- readRDS("expression_matrix.rds")

# Load phenodata (from 3_phenodata_generation.R)
phenodata <- read.csv("phenodata.csv", stringsAsFactors = FALSE)

# Ensure sample order matches columns
stopifnot(identical(colnames(expression_matrix), phenodata$Sample))

# -----------------------------
# 2. Design Matrix for Limma
# -----------------------------
phenodata$Tissue <- factor(phenodata$Tissue, levels = c("DRG", "SN"))

# Create design matrix (model without intercept for each tissue)
design <- model.matrix(~ 0 + Tissue, data = phenodata)
colnames(design) <- c("DRG", "SN")

print("Design matrix:")
print(design)

# -----------------------------
# 3. Linear Modeling with Limma
# -----------------------------
fit <- lmFit(expression_matrix, design)

# Define contrasts (DRG vs SN)
contrast.matrix <- makeContrasts(DRG_vs_SN = DRG - SN, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# -----------------------------
# 4. DEG Extraction Function
# -----------------------------
extract_DEGs <- function(fit2, coef = 1, logFC_cutoff = 1, adjP_cutoff = 0.05,
                         full_results_file = "all_genes_results.csv",
                         deg_file = "DEG_results.csv",
                         volcano_file = "volcano_plot.png") {
  
  # Extract all genes
  results <- topTable(fit2, coef = coef, number = Inf, adjust.method = "BH")
  
  # Save all results
  write.csv(results, full_results_file, row.names = TRUE)
  
  # Filter DEGs
  deg_results <- results %>%
    filter(adj.P.Val < adjP_cutoff & abs(logFC) > logFC_cutoff)
  
  # Save DEGs
  write.csv(deg_results, deg_file, row.names = TRUE)
  
  # Add classification for volcano plot
  results$Significant <- ifelse(results$adj.P.Val < adjP_cutoff &
                                  abs(results$logFC) > logFC_cutoff,
                                "DEG", "Not DEG")
  
  # Volcano plot
  p <- ggplot(results, aes(x = logFC, y = -log10(P.Value), color = Significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("grey", "red")) +
    theme_minimal() +
    labs(title = "Volcano Plot",
         x = "log2 Fold Change",
         y = "-log10(P-value)") +
    theme(legend.title = element_blank())
  
  ggsave(volcano_file, p, width = 6, height = 5, dpi = 300)
  
  # Print summary
  cat("Total genes analyzed:", nrow(results), "\n")
  cat("Significant DEGs (adj.P.Val <", adjP_cutoff,
      "and |logFC| >", logFC_cutoff, "):", nrow(deg_results), "\n")
  cat("Results saved to:", full_results_file, "and", deg_file, "\n")
  cat("Volcano plot saved to:", volcano_file, "\n")
  
  return(deg_results)
}

# -----------------------------
# 5. Extract DEGs and Save Outputs
# -----------------------------
deg_results <- extract_DEGs(fit2, coef = 1)

# Optionally preview top DEGs
head(deg_results)

