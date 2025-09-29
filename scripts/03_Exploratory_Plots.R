if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# List everything you need (Bioc + CRAN)
needs <- c(
  "TCseq", "clusterProfiler", "org.Rn.eg.db", "AnnotationDbi", # Bioc
  "GEOquery", "limma",                                         # Bioc
  "pheatmap", "ggplot2", "tidyverse", "janitor"                # CRAN
)


to_install <- needs[!vapply(needs, requireNamespace, logical(1), quietly = TRUE)]

if (length(to_install)) {
  # BiocManager::install can handle both Bioc and CRAN when repos are set
  BiocManager::install(to_install, update = FALSE, ask = FALSE)
}

library(tidyverse)
library(TCseq)
library(clusterProfiler)
library(org.Rn.eg.db)
library(AnnotationDbi)
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(janitor)
#creae output directory for the figures provided 
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)

#create the annotated table and save it to be used in other steps
get_annotated <- function(fit2, coef_name, annot, out_stub) {
  tt <- limma::topTable(fit2, coef = coef_name, number = Inf, adjust = "fdr") |>
    tibble::rownames_to_column("ProbeID") |>
    dplyr::left_join(annot, by = "ProbeID") |>
    dplyr::mutate(sig = adj.P.Val < 0.05 & abs(logFC) >= 1)
  
  readr::write_csv(tt, file = file.path("outputs/tables", paste0(out_stub, "_DEGs_annotated.csv")))
  saveRDS(tt, file = file.path("outputs/tables", paste0(out_stub, "_DEGs_annotated.rds")))
  tt
}

#create a volacano plot
#log2FC vs -log10 FDR-adjusted P-value
#red denotes a significant DEGs
plot_volcano <- function(tt_annot, title, out_png) {
  p <- ggplot2::ggplot(tt_annot, ggplot2::aes(logFC, -log10(adj.P.Val), color = sig)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.2) +
    ggplot2::scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red"),
                                labels = c("Not sig", "adj.P.Val < 0.05 & |log2FC| ≥ 1"),
                                name = "DEG") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "log2 Fold Change", y = "-log10 adj P-value")
  ggplot2::ggsave(out_png, p, width = 7, height = 5, dpi = 300)
  p
}

#create heatmap of top-N probes
#50 in this case
#add annotations 
plot_heatmap_topN <- function(tt_annot, expr, pheno, N = 50, out_png) {
  probes <- intersect(tt_annot$ProbeID, rownames(expr))
  topn   <- head(probes, N)
  stopifnot(all(colnames(expr) == pheno$Sample))
  mat    <- expr[topn, , drop = FALSE]
  
  ann <- pheno |>
    dplyr::transmute(Sample = Sample,
                     tissue = ifelse(is.na(tissue), "Unknown", tissue),
                     time   = ifelse(is.na(time),   "Unknown", time)) |>
    dplyr::distinct(Sample, .keep_all = TRUE) |>
    tibble::column_to_rownames("Sample")
  
  ann$tissue <- droplevels(factor(ann$tissue))
  ann$time   <- droplevels(factor(ann$time))
  stopifnot(identical(colnames(mat), rownames(ann)))
  
  grDevices::png(out_png, width = 1200, height = 900, res = 150)
  pheatmap::pheatmap(
    mat, scale = "row", annotation_col = ann,
    show_rownames = FALSE, show_colnames = FALSE
  )
  grDevices::dev.off()
}

#create MD (MA) plot from limma for a contrast
plot_md <- function(fit2, coef_name, out_png) {
  grDevices::png(out_png, width = 800, height = 600, res = 120)
  limma::plotMD(fit2, coef = coef_name, status = NULL, main = paste("MD plot:", coef_name))
  abline(h = 0, col = "grey60", lty = 2)
  grDevices::dev.off()
}

#create QQ plot of  t-statistics
plot_qq <- function(fit2, coef_name, out_png) {
  tvals <- fit2$t[, coef_name]
  df    <- fit2$df.total[1]
  grDevices::png(out_png, width = 800, height = 600, res = 120)
  limma::qqt(tvals, df = df, pch = 20, cex = 0.6, main = paste("QQ plot:", coef_name))
  grDevices::dev.off()
}

#create PCA colored by tissue and time 
plot_pca <- function(expr, pheno, out_png) {
  expr_t <- t(expr)
  pc     <- prcomp(expr_t, scale. = TRUE)
  pdat   <- dplyr::tibble(
    Sample = rownames(pc$x),
    PC1 = pc$x[, 1], PC2 = pc$x[, 2]
  ) |>
    dplyr::left_join(pheno, by = "Sample")
  
  p <- ggplot2::ggplot(pdat, ggplot2::aes(PC1, PC2, color = tissue, shape = time)) +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "PCA of Expression", color = "Tissue", shape = "Time")
  
  ggplot2::ggsave(out_png, p, width = 7, height = 5, dpi = 300)
  p
}

#create Sample–sample distance heatmap (another QC view)
plot_sample_distance <- function(expr, pheno, out_png) {
  d <- dist(t(expr))                       # trasnpose the samples
  m <- as.matrix(d)
  # make annotation for columns/rows (same samples and order)
  ann <- pheno |>
    dplyr::transmute(Sample = Sample, tissue = tissue, time = time) |>
    tibble::column_to_rownames("Sample")
  ann <- ann[colnames(m), , drop = FALSE]
  grDevices::png(out_png, width = 1200, height = 1000, res = 150)
  pheatmap::pheatmap(m, clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     annotation_col = ann, annotation_row = ann,
                     show_rownames = FALSE, show_colnames = FALSE,
                     main = "Sample–sample distance (Euclidean)")
  grDevices::dev.off()
}

##########################
# RUN FOR BOTH CONTRASTS
##########################
#below will look at the differnces in 1d and 0d
#this will be for DRG and SN respectively

# DRG 1d vs 0d
drg_tt <- get_annotated(fit2, "DRG_1d_vs_0d", annot, "DRG_1d_vs_0d")
plot_volcano(drg_tt, "Volcano: DRG 1d vs 0d", "outputs/figures/volcano_DRG_1d_vs_0d.png")
plot_heatmap_topN(drg_tt, expr, pheno, 50, "outputs/figures/heatmap_top50_DRG_1d_vs_0d.png")
plot_md(fit2, "DRG_1d_vs_0d", "outputs/figures/MD_DRG_1d_vs_0d.png")
plot_qq(fit2, "DRG_1d_vs_0d", "outputs/figures/QQ_DRG_1d_vs_0d.png")

# SN 1d vs 0d
sn_tt  <- get_annotated(fit2, "SN_1d_vs_0d",  annot, "SN_1d_vs_0d")
plot_volcano(sn_tt,  "Volcano: SN 1d vs 0d",  "outputs/figures/volcano_SN_1d_vs_0d.png")
plot_heatmap_topN(sn_tt, expr, pheno, 50, "outputs/figures/heatmap_top50_SN_1d_vs_0d.png")
plot_md(fit2, "SN_1d_vs_0d", "outputs/figures/MD_SN_1d_vs_0d.png")
plot_qq(fit2, "SN_1d_vs_0d", "outputs/figures/QQ_SN_1d_vs_0d.png")

# global pca to check data quality 
plot_pca(expr, pheno, "outputs/figures/PCA_tissue_time.png")
plot_sample_distance(expr, pheno, "outputs/figures/sample_distance_heatmap.png")

# create Venn diagram overlap of significant DEGs between contrasts
sig_drg <- drg_tt$ProbeID[drg_tt$sig]
sig_sn  <- sn_tt$ProbeID[sn_tt$sig]
venn_mat <- cbind(DRG_1d_vs_0d = rownames(expr) %in% sig_drg,
                  SN_1d_vs_0d  = rownames(expr) %in% sig_sn)
grDevices::png("outputs/figures/venn_DRG_SN.png", width = 800, height = 600, res = 120)
limma::vennDiagram(vennCounts(venn_mat), circle.col = c("tomato", "royalblue"))
grDevices::dev.off()

