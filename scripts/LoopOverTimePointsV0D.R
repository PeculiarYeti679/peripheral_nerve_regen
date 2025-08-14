# ==========================================
# AUTO: all timepoints vs 0d per tissue
# ==========================================

library(limma)
library(tidyverse)
library(pheatmap)
library(janitor)
library(GEOquery)

# ---------- Inputs assumed ----------
# expr  : expression matrix (probes x samples)
# pheno : data.frame with columns Sample, tissue, time
stopifnot(all(colnames(expr) == pheno$Sample))

# ---------- Folders ----------
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/tables",  showWarnings = FALSE, recursive = TRUE)

# ---------- Design (only existing levels) ----------
pheno$tissue <- droplevels(factor(pheno$tissue))
pheno$time   <- droplevels(factor(pheno$time))

design <- model.matrix(~ 0 + tissue:time, data = pheno)
colnames(design) <- make.names(colnames(design))

# ---------- Build contrasts dynamically ----------
times   <- setdiff(levels(pheno$time), "0d")         # all non-baseline timepoints present
tissues <- levels(pheno$tissue)                       # DRG, SN (from your data)

# Expressions like: "DRG_1d_vs_0d = tissueDRG.time1d - tissueDRG.time0d"
contrast_specs <- unlist(lapply(tissues, function(ti) {
  paste0(ti, "_", times, "_vs_0d = tissue", ti, ".time", times, " - tissue", ti, ".time0d")
}))
# Keep only contrasts whose terms actually exist in the design (defensive)
cols_in_design <- colnames(design)
is_valid <- function(expr_str) {
  parts <- strsplit(gsub(" ", "", strsplit(expr_str, "=")[[1]][2]), "-|\\+")[[1]]
  all(parts %in% cols_in_design)
}
contrast_specs <- contrast_specs[vapply(contrast_specs, is_valid, logical(1L))]

stopifnot(length(contrast_specs) > 0)

contrast_matrix <- do.call(makeContrasts, c(list(levels = design), as.list(contrast_specs)))

# ---------- Fit ----------
fit  <- lmFit(expr, design)
fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))

# ---------- Annotation (GPL7294) ----------
gpl      <- GEOquery::getGEO("GPL7294", AnnotGPL = TRUE)
gpl_tbl  <- Table(gpl) |> as_tibble() |> clean_names()
annot <- gpl_tbl %>%
  transmute(
    ProbeID = id,
    SYMBOL  = na_if(gene_symbol, ""),
    ENTREZ  = ifelse(is.na(gene), NA_character_, as.character(gene)), # 'gene' is Entrez ID on this platform
    GENE    = coalesce(gene_name, description)
  ) %>% distinct(ProbeID, .keep_all = TRUE)

# ---------- Helpers ----------
get_annotated <- function(fit2, coef_name, annot, out_stub) {
  tt <- limma::topTable(fit2, coef = coef_name, number = Inf, adjust = "fdr") |>
    tibble::rownames_to_column("ProbeID") |>
    dplyr::left_join(annot, by = "ProbeID") |>
    dplyr::mutate(sig = adj.P.Val < 0.05 & abs(logFC) >= 1)
  readr::write_csv(tt, file.path("outputs/tables", paste0(out_stub, "_DEGs_annotated.csv")))
  saveRDS(tt,        file.path("outputs/tables", paste0(out_stub, "_DEGs_annotated.rds")))
  tt
}

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

# ---------- Run for every contrast ----------
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

# ---------- Optional: global QC once ----------
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

# Sample–sample distance
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

