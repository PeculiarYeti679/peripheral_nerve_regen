#Due to a change in the phenotype data in step 05 there was 
#an update to data that this script uses so step 05 must be ran first. 


library(limma)
library(tidyverse)
library(pheatmap)
library(janitor)
library(GEOquery)


# expr: expression matrix (probes x samples)
# pheno: data.frame with columns Sample, tissue, time
# check to see of the sample names match before continuing
stopifnot(all(colnames(expr) == pheno$Sample))

# output folder
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/tables",  showWarnings = FALSE, recursive = TRUE)

# drop data that does not contain tissue and time
pheno$tissue <- droplevels(factor(pheno$tissue))
pheno$time   <- droplevels(factor(pheno$time))

#create design matrix using tissue and time
design <- model.matrix(~ 0 + tissue:time, data = pheno)
colnames(design) <- make.names(colnames(design))

# create the contrasts at 0d
times   <- setdiff(levels(pheno$time), "0d")         # all non-baseline timepoints present
tissues <- levels(pheno$tissue)                       # DRG, SN (from your data)

# looking for the contrasts where like 0d
contrast_specs <- unlist(lapply(tissues, function(ti) {
  paste0(ti, "_", times, "_vs_0d = tissue", ti, ".time", times, " - tissue", ti, ".time0d")
}))
# filter the contrasts for where there are data
cols_in_design <- colnames(design)
is_valid <- function(expr_str) {
  parts <- strsplit(gsub(" ", "", strsplit(expr_str, "=")[[1]][2]), "-|\\+")[[1]]
  all(parts %in% cols_in_design)
}
contrast_specs <- contrast_specs[vapply(contrast_specs, is_valid, logical(1L))]

stopifnot(length(contrast_specs) > 0)

contrast_matrix <- do.call(makeContrasts, c(list(levels = design), as.list(contrast_specs)))

# create linear model 
fit  <- lmFit(expr, design)
# model with Bayes shrinkage to standard error
# this helps with stability in variance
fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))

# annotate GPL7294
#normalize the column names and allows to from probeID to 
#symbol, entrez and gene name
gpl      <- GEOquery::getGEO("GPL7294", AnnotGPL = TRUE)
gpl_tbl  <- Table(gpl) |> as_tibble() |> clean_names()
annot <- gpl_tbl %>%
  transmute(
    ProbeID = id,
    SYMBOL  = na_if(gene_symbol, ""),
    ENTREZ  = ifelse(is.na(gene), NA_character_, as.character(gene)), # 'gene' is Entrez ID on this platform
    GENE    = coalesce(gene_name, description)
  ) %>% distinct(ProbeID, .keep_all = TRUE)

# get ranked table
# looking for signifcance value dj.P.Val < 0.05 and  abs(logFC) >= 1
get_annotated <- function(fit2, coef_name, annot, out_stub) {
  tt <- limma::topTable(fit2, coef = coef_name, number = Inf, adjust = "fdr") |>
    tibble::rownames_to_column("ProbeID") |>
    dplyr::left_join(annot, by = "ProbeID") |>
    dplyr::mutate(sig = adj.P.Val < 0.05 & abs(logFC) >= 1)
  readr::write_csv(tt, file.path("outputs/tables", paste0(out_stub, "_DEGs_annotated.csv")))
  saveRDS(tt,        file.path("outputs/tables", paste0(out_stub, "_DEGs_annotated.rds")))
  tt
}

#far right or left indicate a significant change
plot_volcano <- function(tt_annot, title, out_png) {
  p <- ggplot(tt_annot, aes(logFC, -log10(adj.P.Val), color = sig)) +
    geom_point(alpha = 0.7, size = 1.2) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red"),
                       labels = c("Not sig", "adj.P.Val < 0.05 & |log2FC| ≥ 1"),
                       name = "DEG") +
    theme_minimal() +
    labs(title = title, x = "log2 Fold Change", y = "-log10 adj P-value")
  ggsave(out_png, p, width = 7, height = 5, dpi = 300)
}

#create heatmap to see if clusters of samples group together
plot_heatmap_topN <- function(tt_annot, expr, pheno, N, out_png) {
  probes <- intersect(tt_annot$ProbeID, rownames(expr))
  topn   <- head(probes, N)
  if (length(topn) < 2) return(invisible(NULL))
  stopifnot(all(colnames(expr) == pheno$Sample))
  mat <- expr[topn, , drop = FALSE]
  
  ann <- pheno |>
    transmute(Sample = Sample,
              tissue = ifelse(is.na(tissue), "Unknown", tissue),
              time   = ifelse(is.na(time),   "Unknown", time)) |>
    distinct(Sample, .keep_all = TRUE) |>
    column_to_rownames("Sample")
  ann$tissue <- droplevels(factor(ann$tissue))
  ann$time   <- droplevels(factor(ann$time))
  stopifnot(identical(colnames(mat), rownames(ann)))
  
  png(out_png, width = 1200, height = 900, res = 150)
  pheatmap(mat, scale = "row", annotation_col = ann,
           show_rownames = FALSE, show_colnames = FALSE)
  dev.off()
}

plot_md <- function(fit2, coef_name, out_png) {
  png(out_png, width = 800, height = 600, res = 120)
  limma::plotMD(fit2, coef = coef_name, main = paste("MD plot:", coef_name))
  abline(h = 0, col = "grey60", lty = 2)
  dev.off()
}

plot_qq <- function(fit2, coef_name, out_png) {
  tvals <- fit2$t[, coef_name]
  df    <- fit2$df.total[1]
  png(out_png, width = 800, height = 600, res = 120)
  limma::qqt(tvals, df = df, pch = 20, cex = 0.6, main = paste("QQ plot:", coef_name))
  dev.off()
}

# running these for all contrasts
all_contrasts <- colnames(contrast_matrix)

for (cn in all_contrasts) {
  # out_stub like "DRG_1d_vs_0d"
  out_stub <- cn
  
  # 1) annotate + save tables
  tt <- get_annotated(fit2, cn, annot, out_stub)
  
  # 2) figures
  plot_volcano(tt, paste("Volcano:", cn), file.path("outputs/figures", paste0("volcano_", out_stub, ".png")))
  plot_heatmap_topN(tt, expr, pheno, 50,  file.path("outputs/figures", paste0("heatmap_top50_", out_stub, ".png")))
  plot_md(fit2, cn, file.path("outputs/figures", paste0("MD_", out_stub, ".png")))
  plot_qq(fit2, cn, file.path("outputs/figures", paste0("QQ_", out_stub, ".png")))
}

# QC data 
# PCA
expr_t <- t(expr)
pc     <- prcomp(expr_t, scale. = TRUE)
pdat   <- tibble(Sample = rownames(pc$x), PC1 = pc$x[,1], PC2 = pc$x[,2]) |>
  left_join(pheno, by = "Sample")
p <- ggplot(pdat, aes(PC1, PC2, color = tissue, shape = time)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_minimal() +
  labs(title = "PCA of Expression", color = "Tissue", shape = "Time")
ggsave("outputs/figures/PCA_tissue_time.png", p, width = 7, height = 5, dpi = 300)

# sample–sample distance
d <- dist(t(expr)); m <- as.matrix(d)
annS <- pheno |>
  transmute(Sample = Sample, tissue = tissue, time = time) |>
  column_to_rownames("Sample")
annS <- annS[colnames(m), , drop = FALSE]
png("outputs/figures/sample_distance_heatmap.png", width = 1200, height = 1000, res = 150)
pheatmap(m, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         annotation_col = annS, annotation_row = annS,
         show_rownames = FALSE, show_colnames = FALSE,
         main = "Sample–sample distance (Euclidean)")
dev.off()



library(dplyr)
library(readr)

# collect GO terms across clusters
go_files <- list.files("final", pattern="GO_BP_cluster_\\d+\\.csv", full.names=TRUE)
go_all <- bind_rows(lapply(go_files, read_csv, show_col_types = FALSE), .id="cluster_id")
colnames(go_all)
head(go_all)

# add cluster number
go_all <- go_all %>%
  mutate(Cluster = as.integer(gsub("\\D","", cluster_id)))

# top 5 terms per cluster
go_top5 <- go_all %>%
  group_by(Cluster) %>%
  arrange(p.adjust) %>%
  slice_head(n=5) %>%
  ungroup()

write_csv(go_top5, "final/GO_BP_top5_summary.csv")


# collect KEGG pathways
kegg_files <- list.files("final", pattern="KEGG_cluster_\\d+\\.csv", full.names=TRUE)
kegg_all <- bind_rows(lapply(kegg_files, read_csv, show_col_types = FALSE), .id="cluster_id")

kegg_all <- kegg_all %>%
  mutate(Cluster = as.integer(gsub("\\D","", cluster_id)))

kegg_top5 <- kegg_all %>%
  group_by(Cluster) %>%
  arrange(p.adjust) %>%
  slice_head(n=5) %>%
  ungroup()

write_csv(kegg_top5, "final/KEGG_top5_summary.csv")
library(dplyr)
library(stringr)
library(readr)


go_all <- read_csv("final/GO_BP_top5_summary.csv", show_col_types = FALSE)

# define theme keywords
themes <- list(
  Immune        = c("immune", "cytokine", "leukocyte", "macrophage", "inflammatory", "antigen", "defense"),
  CellCycle     = c("cell cycle", "mitotic", "replication", "checkpoint", "division", "proliferation"),
  AxonRepair    = c("axon", "synapse", "neuro", "dendrite", "neuron", "synaptic", "plasticity", "guidance"),
  Metabolism    = c("metabolic", "oxidative", "mitochond", "respiratory", "glycolysis", "lipid", "glucose"),
  Stress        = c("stress", "response", "MAPK", "apoptosis", "DNA damage", "hypoxia")
)

# function to assign theme
assign_theme <- function(term) {
  term_l <- tolower(term)
  for (th in names(themes)) {
    if (any(str_detect(term_l, themes[[th]]))) return(th)
  }
  return("Other")
}

# apply to GO terms
go_all <- go_all %>%
  mutate(Theme = sapply(Description, assign_theme))

# save
write_csv(go_all, "final/GO_BP_top5_withThemes.csv")
list.files("final", pattern="GO_BP_cluster_\\d+\\.csv")

# same for KEGG
kegg_all <- read_csv("final/KEGG_top5_summary.csv", show_col_types = FALSE) %>%
  mutate(Theme = sapply(Description, assign_theme))
write_csv(kegg_all, "final/KEGG_top5_withThemes.csv")
