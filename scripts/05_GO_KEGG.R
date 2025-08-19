if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.21", ask = FALSE)

library(tidyverse)
library(TCseq)
library(clusterProfiler)
library(org.Rn.eg.db)
library(pheatmap)
library(ggplot2)

#create directory for final findings
dir.create("final", showWarnings = FALSE)

#load in data
expression_matrix <- readRDS("data/processed/expression_matrix_processed.rds")
phenodata <- read.csv("data/meta/phenodata.csv", stringsAsFactors = FALSE)

#Check if data is of the right type
stopifnot(is.matrix(expression_matrix) || is.data.frame(expression_matrix))
stopifnot(is.data.frame(phenodata))


#standardize the data
expression_matrix <- as.matrix(expression_matrix)

colnames(phenodata)

#verify that columns are aligned
all(colnames(expression_matrix) %in% phenodata$Sample)
sum(!colnames(expression_matrix) %in% phenodata$Sample)
phenodata %>% head() %>% as_tibble()

#could use DRG
tissue_of_interest <- "SN"

#get all SN from the phenodata
pd_tissue <- phenodata %>% filter(tissue == tissue_of_interest)

#times in the data
time_levels <-c("0d","1d","4d","7d","14d")

#force data to be in correct order
pd_tissue$time <- factor(pd_tissue$time,
                         levels = intersect(time_levels, unique(pd_tissue$time)),
                         ordered = TRUE)
pd_tissue <- pd_tissue %>% arrange(time)
expr_tissue <- expression_matrix[, pd_tissue$Sample, drop = FALSE]


#verify data
table(pd_tissue$Time)
dim(expr_tissue)

#remove NA data
expr_tissue <- expr_tissue[rowSums(is.na(expr_tissue)) == 0, , drop=FALSE]

var_genes <- apply(expr_tissue, 1, var)
var_genes <- var_genes[is.finite(var_genes) & !is.na(var_genes)]

#keep 3000 genes
keep_n <- min(3000, length(var_genes))
sig_genes <- names(sort(var_genes, decreasing = TRUE))[1:keep_n]

expr_deg <- expr_tissue[rownames(expr_tissue) %in% sig_genes, , drop = FALSE]
dim(expr_deg)


set.seed(42)
k_clusters <- 6 

#km = k means
tc <- timeclust(expr_deg, algo = "km", k = k_clusters, standardize = TRUE)

#get cluster labels
clusters <- tryCatch(tc@cluster, error = function(e) NULL)
#check if clusters is null and get from the time clusters
if (is.null(clusters)) clusters <- tc$cluster 

#create cluster dataframe
cluster_df <- tibble(Gene = rownames(expr_deg), Cluster = as.integer(clusters))
write.csv(cluster_df, "final/tcseq_clusters.csv", row.names = FALSE)

#check the clusters
table(cluster_df$Cluster)

k_ids <- sort(unique(cluster_df$Cluster))
timeclustplot(tc, value="z-score", cols=3,)
n_panels <- length(k_ids)
n_cols <- min(3, n_panels)

# create the trajectory plot 
png("final/tcseq_cluster_trajectories.png", width = 1200, height = 800, res = 150)
#this timecluster is hard to read and due to the amount of data
timeclustplot(tc, value = "z-score", cols = n_cols)
dev.off()

#using average z-score per cluster per timepoint
#this should make a an easier to interrupt tracjectory chart

# cluster_df built earlier: Gene + Cluster
expr_z <- t(scale(t(expr_deg)))

# compute per-cluster centroid
centroids <- sapply(sort(unique(cluster_df$Cluster)), function(k) {
  rows <- which(cluster_df$Cluster == k)
  if (length(rows) == 0) rep(NA_real_, ncol(expr_z)) else colMeans(expr_z[rows, , drop = FALSE])
})
colnames(centroids) <- paste0("C", sort(unique(cluster_df$Cluster)))

cent_df <- as.data.frame(centroids) |>
  dplyr::mutate(Sample = colnames(expr_z)) |>
  dplyr::left_join(pd_tissue |> dplyr::select(Sample, time, time_hours), by = "Sample") |>
  tidyr::pivot_longer(-c(Sample, time, time_hours),
                      names_to = "Cluster", values_to = "Z") |>
  dplyr::group_by(Cluster, time, time_hours) |>
  dplyr::summarize(Z = mean(Z, na.rm=TRUE), .groups="drop") |>
  dplyr::arrange(time_hours)

ggplot(cent_df, aes(time_hours, Z, group=Cluster, color=Cluster)) +
  geom_line(size=1.2) + geom_point(size=2) +
  facet_wrap(~ Cluster, scales="free_y") +
  scale_x_continuous(breaks = unique(cent_df$time_hours),
                     labels = unique(cent_df$time)) +
  labs(title=paste("Cluster trajectories (", tissue_of_interest, ")", sep=""),
       x="Time", y="Mean z-score") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("final/tcseq_cluster_trajectories_faceted.png", width=10, height=7, dpi=300)





# create heatmaps ordered by cluters
expr_z <- t(scale(t(expr_deg)))
ord_idx <- order(cluster_df$Cluster)
#get rownames and align them
mat <- expr_z[ord_idx, , drop = FALSE]
ann_row <- data.frame(Cluster = factor(cluster_df$Cluster[ord_idx]))
rownames(ann_row) <- rownames(mat) 

pheatmap::pheatmap(
  mat,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE,
  annotation_row = ann_row,
  filename = "final/tcseq_heatmap.png",
  width = 9, height = 10
)

# Faceted cluster-mean trajectories (cleaner captioning)
centroids <- sapply(1:k_clusters, function(k) {
  rows <- which(cluster_df$Cluster == k)
  if (length(rows) == 0) rep(NA_real_, ncol(expr_z)) else colMeans(expr_z[rows, , drop = FALSE])
})
colnames(centroids) <- paste0("C", 1:k_clusters)

cent_df <- as.data.frame(centroids) %>%
  mutate(Sample = colnames(expr_z)) %>%
  left_join(pd_tissue %>% dplyr::select(Sample, time), by = "Sample") %>%
  pivot_longer(-c(Sample, time), names_to = "Cluster", values_to = "Z") %>%
  group_by(Cluster, time) %>%
  summarize(Z = mean(Z, na.rm=TRUE), .groups="drop") %>%
  mutate(Time = factor(time, levels = levels(pd_tissue$time), ordered = TRUE))

ggplot(cent_df, aes(Time, Z, group = Cluster)) +
  geom_line() + geom_point() +
  facet_wrap(~ Cluster, scales = "free_y") +
  labs(title = paste("TCseq cluster trajectories -", tissue_of_interest),
       x = "Time", y = "Mean z-score") +
  theme_bw()
ggsave("final/tcseq_cluster_trajectories_faceted.png", width = 10, height = 7, dpi = 300)

##########

#GO/KEGG Enrichment 
library(GEOquery)
library(limma)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(clusterProfiler)
library(org.Rn.eg.db)   # Rattus norvegicus
library(AnnotationDbi)

#helper functions to get GPL COLUMNS
find_first <- function(pattern, cn) {
  hit <- grep(pattern, cn, ignore.case = TRUE, value = TRUE)
  if (length(hit)) hit[1] else NA_character_
}
split_multi <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\\s*///\\s*|\\s*[,;|]\\s*|\\s+/\\s+", ";")
  ifelse(nzchar(x), x, NA_character_)
}
normalize_ens <- function(x) {
  x <- toupper(as.character(x))
  sub("\\.\\d+$", "", x)  # strip version suffix .1, .2, ...
}
gse <- getGEO("GSE30165", GSEMatrix = TRUE)
if (is.list(gse)) {
  esets <- gse
} else {
  esets <- list(gse)
}

gpl_id <- annotation(eset)
gpl <- getGEO(gpl_id)
ann <- Table(gpl)
cn <- colnames(ann)

probe_col  <- if (!is.na(find_first("^ID$", cn))) "ID" else if (!is.na(find_first("^ID_REF$", cn))) "ID_REF" else cn[1]
symbol_col <- find_first("^gene.*symbol$|(^|_)symbol($|_)|symbol", cn)
ens_col    <- find_first("ensembl", cn)
entrez_col <- find_first("entrez.*(gene)?\\s*id|geneid|entrezid|entrez_id", cn)

cat("Columns picked:\n  Probe:", probe_col, "\n  SYMBOL:", symbol_col, "\n  ENSEMBL:", ens_col, "\n  ENTREZ:", entrez_col, "\n")


annot_raw <- ann |>
  transmute(
    ProbeID   = .data[[probe_col]],
    SYMBOL    = if (!is.na(symbol_col)) split_multi(.data[[symbol_col]]) else NA,
    ENSEMBL   = if (!is.na(ens_col))    split_multi(.data[[ens_col]])    else NA,
    ENTREZID0 = if (!is.na(entrez_col)) split_multi(.data[[entrez_col]]) else NA
  ) |>
  mutate(
    SYMBOL    = ifelse(!is.na(SYMBOL),  str_split_fixed(SYMBOL,  ";", 2)[,1], NA),
    ENSEMBL   = ifelse(!is.na(ENSEMBL), normalize_ens(str_split_fixed(ENSEMBL, ";", 2)[,1]), NA),
    ENTREZID0 = ifelse(!is.na(ENTREZID0), str_split_fixed(ENTREZID0, ";", 2)[,1], NA)
  ) |>
  distinct(ProbeID, .keep_all = TRUE)

# clean empties
for (v in c("SYMBOL","ENSEMBL","ENTREZID0")) {
  annot_raw[[v]][annot_raw[[v]] %in% c("", "NA", "NULL", "null")] <- NA
}

has_sym <- !is.na(annot_raw$SYMBOL)
entrez_from_symbol <- rep(NA_character_, nrow(annot_raw))
if (any(has_sym)) {
  m <- AnnotationDbi::mapIds(org.Rn.eg.db,
                             keys    = annot_raw$SYMBOL[has_sym],
                             keytype = "SYMBOL",
                             column  = "ENTREZID",
                             multiVals = "first")
  entrez_from_symbol[has_sym] <- unname(m[annot_raw$SYMBOL[has_sym]])
}

is_rat_gene <- !is.na(annot_raw$ENSEMBL) & grepl("^ENSRNOG", annot_raw$ENSEMBL)
entrez_from_ens <- rep(NA_character_, nrow(annot_raw))
if (any(is_rat_gene)) {
  ens_keys <- unique(annot_raw$ENSEMBL[is_rat_gene])
  m2 <- AnnotationDbi::mapIds(org.Rn.eg.db,
                              keys    = ens_keys,
                              keytype = "ENSEMBL",
                              column  = "ENTREZID",
                              multiVals = "first")
  idx <- which(is_rat_gene)
  ens_vals <- annot_raw$ENSEMBL[idx]
  entrez_from_ens[idx] <- unname(m2[ens_vals])
}

annot <- annot_raw |>
  mutate(ENTREZID = coalesce(ENTREZID0, entrez_from_symbol, entrez_from_ens)) |>
  dplyr::select(ProbeID, SYMBOL, ENSEMBL, ENTREZID)

# Coverage summary
cov <- c(
  platform_entrez = sum(!is.na(annot_raw$ENTREZID0)),
  via_SYMBOL      = sum(!is.na(entrez_from_symbol)),
  via_ENSEMBL     = sum(!is.na(entrez_from_ens)),
  final_nonNA     = sum(!is.na(annot$ENTREZID)),
  probes_total    = nrow(annot)
)
print(cov)


#phenodata 
pd <- pData(eset) |>
  as.data.frame() |>
  tibble::rownames_to_column("Sample") |>
  tibble::as_tibble()

# columns like characteristics_ch1, characteristics_ch1.1, etc.
char_cols <- grep("^characteristics.*ch1", names(pd), ignore.case = TRUE, value = TRUE)

pd_long <- pd |>
  dplyr::select(
    Sample,
    title,
    source_name_ch1 = dplyr::any_of("source_name_ch1"),
    dplyr::all_of(char_cols)
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::starts_with("characteristics"),
    names_to = "char",
    values_to = "text",
    values_drop_na = TRUE
  ) |>
  # split on the first ":"; keep the rest in value
  tidyr::separate(
    col   = text,
    into  = c("key","value"),
    sep   = ":[[:space:]]*",
    extra = "merge",
    fill  = "right"
  ) |>
  dplyr::mutate(
    key   = stringr::str_to_lower(stringr::str_trim(key)),
    value = stringr::str_trim(value)
  ) |>
  dplyr::filter(!is.na(key), !is.na(value)) |>
  dplyr::distinct(Sample, key, .keep_all = TRUE)

pd_wide <- pd_long |>
  tidyr::pivot_wider(
    id_cols    = Sample,
    names_from = key,
    values_from= value,
    values_fn  = \(x) x[1]
  )
# Heuristics to derive tissue, time, status from any available columns
all_pd <- pd |>
  left_join(pd_wide, by = "Sample") |>
  mutate(
    tissue = coalesce(
      .data[["tissue"]], .data[["tissue type"]], .data[["source_name_ch1"]],
      ifelse(str_detect(title, regex("DRG|dorsal root", ignore_case = TRUE)), "DRG",
             ifelse(str_detect(title, regex("sciatic|SN", ignore_case = TRUE)), "SN", NA))
    ),
    time_raw = coalesce(.data[["time"]], .data[["time point"]], .data[["timepoint"]], .data[["time_point"]], title),
    status = coalesce(.data[["status"]], .data[["treatment"]], .data[["group"]],
                      ifelse(str_detect(title, regex("injur|axotom|lesion", ignore_case = TRUE)), "injured",
                             ifelse(str_detect(title, regex("control|naive|sham|uninjur", ignore_case = TRUE)), "control", NA)))
  )

# Standardize time to tokens like "0d","1d","4d","1h"...
parse_time <- function(x) {
  x <- tolower(x)
  x <- str_replace_all(x, "days?", "d")
  x <- str_replace_all(x, "hours?|hrs?", "h")
  m <- str_match(x, "\\b(\\d+(?:\\.\\d+)?)\\s*(d|h)\\b")
  out <- ifelse(!is.na(m[,1]), paste0(m[,2], m[,3]), NA)
  out
}

all_pd <- all_pd |>
  mutate(
    tissue = toupper(str_replace_all(tissue, c("dorsal.*ganglia"="DRG", "sciatic.*"="SN", "sn"="SN"))),
    time   = parse_time(time_raw),
    status = tolower(status),
    status = case_when(
      str_detect(status, "injur|axotom|lesion") ~ "injured",
      str_detect(status, "control|naive|sham|uninjur") ~ "control",
      TRUE ~ status
    )
  )

# Keep only rows with essentials
pheno <- all_pd |>
  select(Sample, tissue, time, status, title) |>
  distinct()

# Basic sanity
cat("Counts by tissue/time/status:\n")
print(pheno |>
        count(tissue, time, status) |>
        arrange(tissue, time, desc(n)))
stopifnot(all(colnames(exprs(eset)) %in% pheno$Sample))
# reorder pheno to match expression columns
pheno <- pheno[match(colnames(exprs(eset)), pheno$Sample), ]
stopifnot(all(pheno$Sample == colnames(exprs(eset))))