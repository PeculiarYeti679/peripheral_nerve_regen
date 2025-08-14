# ===============================
# Master Pipeline for GSE30165
# Author: (auto-generated)
# Purpose: End-to-end pipeline orchestrating preprocessing, QC, phenodata,
#          PCA, differential expression, timecourse analysis, enrichment,
#          and report export. All outputs go to final/[...].
# ===============================

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(readr)
  library(ggplot2)
  library(pheatmap)
})

# -------------
# Config
# -------------
options(stringsAsFactors = FALSE)
SEED <- 42
set.seed(SEED)

# Input options:
# - If you already downloaded GEO supplementary files, set RAW_DIR to the folder
#   containing the .txt.gz expression files you currently use.
# - Otherwise, set DOWNLOAD_GEO <- TRUE to fetch with GEOquery.
DOWNLOAD_GEO <- FALSE
GSE_ACC <- "GSE30165"
RAW_DIR <- "GSE30165_raw"        # folder containing untarred .txt.gz
PLATFORM <- "Agilent"            # used for labels only

# Output root (created if missing)
OUT <- "final"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dirs <- file.path(OUT, c("data","qc","phenodata","pca","deg","timecourse","enrichment","report"))
invisible(lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE))

# -------------
# Helpers
# -------------
message2 <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), sprintf(...), "\n")

save_plot <- function(p, path, width = 8, height = 6, dpi = 300) {
  ggsave(filename = path, plot = p, width = width, height = height, dpi = dpi, limitsize = FALSE)
}

# -------------
# (1) Preprocessing
# -------------
message2("STEP 1/7: Preprocessing")

if (DOWNLOAD_GEO) {
  suppressPackageStartupMessages(library(GEOquery))
  message2("Downloading GEO supplementary files for %s ...", GSE_ACC)
  getGEOSuppFiles(GSE_ACC, baseDir = RAW_DIR)
  # Note: user may need to untar manually depending on platform specifics.
}

# Load all .txt.gz into a single expression matrix
# Expectation: columns = samples, rows = probes/genes (after merge)
txt_files <- list.files(RAW_DIR, pattern = "\\.txt(\\.gz)?$", full.names = TRUE, recursive = TRUE)
stopifnot(length(txt_files) > 0)

message2("Found %d raw files", length(txt_files))

# Example reader for Agilent-like single channel TXT (user may adapt):
# We will attempt to read each file, keep a common 'ProbeName' or first column as row id.
read_one <- function(fp) {
  df <- suppressWarnings(suppressMessages(read_tsv(fp, guess_max = 100000)))
  # Heuristics to find signal column
  cand_cols <- c("gProcessedSignal","ProcessedSignal","Signal","VALUE","Intensity","gMedianSignal")
  sig_col <- cand_cols[cand_cols %in% names(df)]
  if (length(sig_col) == 0) {
    stop(sprintf("No signal column found in %s. Please adjust cand_cols.", basename(fp)))
  }
  # Heuristics to find probe id
  id_cols <- c("ProbeName","ProbeID","ID","FeatureNum","SystematicName")
  id_col <- id_cols[id_cols %in% names(df)]
  if (length(id_col) == 0) {
    # fallback to row number
    df$ProbeID <- seq_len(nrow(df))
    id_col <- "ProbeID"
  } else {
    id_col <- id_col[1]
  }
  out <- df[, c(id_col, sig_col[1])]
  names(out) <- c("ProbeID", tools::file_path_sans_ext(basename(fp)))
  out
}

lst <- lapply(txt_files, read_one)
expr <- Reduce(function(x,y) full_join(x, y, by = "ProbeID"), lst)
stopifnot("ProbeID" %in% names(expr))

# Make rownames and numeric matrix
rn <- expr$ProbeID
mat <- as.matrix(expr[ , setdiff(names(expr), "ProbeID")])
mode(mat) <- "numeric"
rownames(mat) <- rn

# Basic integrity checks
message2("Expression matrix dims: %d rows x %d cols", nrow(mat), ncol(mat))
stopifnot(!anyNA(mat))

# Save raw matrix
saveRDS(mat, file.path(OUT, "data", "expression_matrix_raw.rds"))

# Log2 transform and quantile normalize with limma
message2("Applying log2 transform and between-array normalization")
mat_log <- log2(mat + 1)
mat_norm <- limma::normalizeBetweenArrays(mat_log, method = "quantile")

saveRDS(mat_norm, file.path(OUT, "data", "expression_matrix_norm.rds"))
write.csv(mat_norm, file.path(OUT, "data", "expression_matrix_norm.csv"), quote = FALSE)

# -------------
# (2) QC
# -------------
message2("STEP 2/7: QC & visualization")

# Boxplot (distribution)
df_long <- as.data.frame(mat_norm) %>%
  rownames_to_column("ProbeID") %>%
  pivot_longer(-ProbeID, names_to = "Sample", values_to = "Expression")

p_box <- ggplot(df_long, aes(x = Sample, y = Expression)) +
  geom_boxplot(outlier.size = 0.2) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = sprintf("%s normalized expression (log2)", PLATFORM))
save_plot(p_box, file.path(OUT, "qc", "boxplot_normalized.png"), width = 12, height = 6)

# Density plot
p_dens <- ggplot(df_long, aes(Expression, group = Sample)) +
  geom_density(alpha = 0.3) +
  theme_bw() +
  labs(title = "Density of normalized expression")
save_plot(p_dens, file.path(OUT, "qc", "density_normalized.png"))

# Sample-sample correlation heatmap
cors <- cor(mat_norm, method = "spearman", use = "pairwise.complete.obs")
png(file.path(OUT, "qc", "sample_correlation_heatmap.png"), width = 1200, height = 1000, res = 150)
pheatmap(cors, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
         main = "Sample Correlation (Spearman)")
dev.off()

# -------------
# (3) Phenodata generation
# -------------
message2("STEP 3/7: Phenodata generation")

samples <- colnames(mat_norm)

# Parse patterns like "DRG_sciatic nerve resection_0d" from sample names if present.
# Users may need to adapt this regex to match actual file naming.
parsed <- tibble(Sample = samples) %>%
  mutate(raw = Sample,
         tissue = str_extract(raw, "^(DRG|SN|SciaticNerve|DorsalRootGanglia)"),
         timepoint = str_extract(raw, "(0\\.5h|1h|3h|6h|9h|1d|4d|7d|14d|0d|3d)"),
         condition = case_when(
           str_detect(raw, "resection|injur|cut|lesion") ~ "Injured",
           TRUE ~ "Control"
         ))

# Fallbacks if regex failed
parsed <- parsed %>%
  mutate(tissue = ifelse(is.na(tissue), "Unknown", tissue),
         timepoint = ifelse(is.na(timepoint), "Unknown", timepoint))

write.csv(parsed, file.path(OUT, "phenodata", "phenodata.csv"), row.names = FALSE)

# -------------
# (4) PCA
# -------------
message2("STEP 4/7: PCA & clustering")

# center-scale by gene
mat_scaled <- t(scale(t(mat_norm)))
pcs <- prcomp(t(mat_scaled), center = TRUE, scale. = FALSE)

var_expl <- pcs$sdev^2 / sum(pcs$sdev^2) * 100
pc_df <- data.frame(Sample = rownames(pcs$x),
                    PC1 = pcs$x[,1],
                    PC2 = pcs$x[,2]) %>%
  left_join(parsed, by = "Sample")

p_pca <- ggplot(pc_df, aes(PC1, PC2, label = Sample, color = tissue, shape = condition)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = sprintf("PCA of %s (PC1 %.1f%%, PC2 %.1f%%)", PLATFORM, var_expl[1], var_expl[2]))
save_plot(p_pca, file.path(OUT, "pca", "pca_scatter.png"))

# Dendrogram / clustering heatmap on top variable genes
vars <- apply(mat_norm, 1, var, na.rm = TRUE)
top_idx <- order(vars, decreasing = TRUE)[seq_len(min(2000, length(vars)))]
png(file.path(OUT, "pca", "heatmap_topvar2000.png"), width = 1200, height = 1200, res = 150)
pheatmap(mat_norm[top_idx, ], show_rownames = FALSE, show_colnames = FALSE,
         main = "Top variable genes (up to 2000)")
dev.off()

# -------------
# (5) Differential Expression (limma)
# -------------
message2("STEP 5/7: Differential expression with limma")

# Build a minimal design using condition (and optionally tissue/timepoint if available)
design_df <- parsed %>%
  mutate(condition = factor(condition),
         tissue = factor(tissue),
         timepoint = factor(timepoint))

# Example design: ~ 0 + condition
# Adjust to: ~ 0 + tissue + condition, or include timepoint for richer model if balanced
design <- model.matrix(~ 0 + condition, data = design_df)
colnames(design) <- levels(design_df$condition)

# Fit with limma
fit <- lmFit(mat_norm, design)

# Example contrast: Injured vs Control (if both exist)
if (all(c("Control","Injured") %in% colnames(design))) {
  contrast.matrix <- makeContrasts(Injured_vs_Control = Injured - Control, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  tt <- topTable(fit2, coef = "Injured_vs_Control", number = Inf, sort.by = "P")
  tt$ProbeID <- rownames(tt)
  write.csv(tt, file.path(OUT, "deg", "DEG_Injured_vs_Control.csv"), row.names = FALSE)

  # Volcano
  p_volcano <- ggplot(tt, aes(logFC, -log10(adj.P.Val))) +
    geom_point(alpha = 0.5) +
    theme_bw() +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(title = "Volcano: Injured vs Control")
  save_plot(p_volcano, file.path(OUT, "deg", "volcano_Injured_vs_Control.png"))
}

# Save full design and fit objects
saveRDS(list(design = design, fit = fit), file.path(OUT, "deg", "limma_fit.rds"))

# -------------
# (6) Timecourse / Trend analysis (basic)
# -------------
message2("STEP 6/7: Timecourse (basic aggregation)")

# A simple approach: compute mean per (tissue, timepoint) and save for downstream clustering
agg_keys <- c("tissue","timepoint")
stopifnot(all(agg_keys %in% names(parsed)))

sample_map <- parsed %>% select(Sample, all_of(agg_keys))

# average expression per group
group_means <- lapply(split(sample_map, interaction(sample_map$tissue, sample_map$timepoint, drop = TRUE)), function(df) {
  s <- df$Sample
  if (length(s) == 1) {
    mat_norm[, s, drop = FALSE]
  } else {
    rowMeans(mat_norm[, s, drop = FALSE])
  }
})
group_means_mat <- do.call(cbind, group_means)
colnames(group_means_mat) <- names(group_means)

saveRDS(group_means_mat, file.path(OUT, "timecourse", "group_means.rds"))
write.csv(group_means_mat, file.path(OUT, "timecourse", "group_means.csv"))

# Optional: soft clustering or maSigPro could be added here.

# -------------
# (7) Enrichment stubs (optional)
# -------------
message2("STEP 7/7: Enrichment (stub, optional)")
# Users can uncomment to run if packages are installed and organism mapping is defined
# suppressPackageStartupMessages({
#   library(clusterProfiler); library(org.Rn.eg.db)
# })
# Example:
# sig <- subset(tt, adj.P.Val < 0.05 & abs(logFC) >= 1)
# genes <- unique(na.omit(sig$ENTREZID))   # requires annotation step
# ego <- enrichGO(gene = genes, OrgDb = org.Rn.eg.db, keyType = "ENTREZID", ont = "BP")
# write.csv(as.data.frame(ego), file.path(OUT, "enrichment", "GO_BP.csv"), row.names = FALSE)

message2("Pipeline complete. Outputs written to 'final/'.")
