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

#get entrez id from the GENE_SYMBOL column of the ann table
entrez_id <- mapIds(org.Rn.eg.db, keys = ann$GENE_SYMBOL,
                   column = "ENTREZID", keytype = "SYMBOL")
#remove any NA 
entrez_id <- entrez_id[!is.na(entrez_id)]

sym_candidates <- colnames(ann)[grepl("symbol", colnames(ann), ignore.case = TRUE)]
sym_col <- if (length(sym_candidates)) sym_candidates[1] else NA_character_
sym_col

probe2sym <- ann %>%
  transmute(
    ProbeID = ID,
    SYMBOL_RAW = .data[[sym_col]]
  ) %>%
  filter(!is.na(SYMBOL_RAW), SYMBOL_RAW != "") %>%
  mutate(SYMBOL = str_split(SYMBOL_RAW, "\\s*///\\s*|\\s*[;,]\\s*")) %>%
  unnest(SYMBOL) %>%
  mutate(SYMBOL = str_trim(SYMBOL)) %>%
  filter(SYMBOL != "") %>%
  distinct(ProbeID, SYMBOL)
sum(rownames(expression_matrix) %in% probe2sym$ProbeID)

head(probe2sym)


expr_gene_mat <- as.data.frame(expression_matrix) %>%
  tibble::rownames_to_column("ProbeID") %>%
  inner_join(probe2sym, by = "ProbeID") %>%
  dplyr::select(SYMBOL, all_of(colnames(expression_matrix))) %>%
  group_by(SYMBOL) %>%
  summarize(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
  tibble::column_to_rownames("SYMBOL") %>%
  as.matrix()

expr_tissue <- expr_gene_mat[, pd_tissue$Sample, drop = FALSE]


expr_tissue <- expr_tissue[rowSums(is.na(expr_tissue)) == 0, , drop = FALSE]
var_genes <- apply(expr_tissue, 1, var)
sig_genes <- names(sort(var_genes[is.finite(var_genes)], decreasing = TRUE))[1:min(3000, length(var_genes))]
expr_deg   <- expr_tissue[rownames(expr_tissue) %in% sig_genes, , drop = FALSE]
dim(expr_deg)   
head(rownames(expr_deg))


set.seed(42)
expr_z <- t(scale(t(expr_deg)))
k_clusters <- 6
km <- kmeans(expr_z, centers = k_clusters, nstart = 25)

cluster_df <- tibble::tibble(
  Gene    = rownames(expr_deg),         # <-- SYMBOLs now
  Cluster = as.integer(km$cluster)
)

centroids <- sapply(1:k_clusters, function(k) colMeans(expr_z[km$cluster == k, , drop = FALSE]))
colnames(centroids) <- paste0("C", 1:k_clusters)

table(cluster_df$Cluster)

library(clusterProfiler); library(org.Rn.eg.db); library(AnnotationDbi)
dir.create("final", showWarnings = FALSE)

# robust SYMBOL -> ENTREZ mapper
symbol_to_entrez <- function(genes) {
  m1 <- AnnotationDbi::mapIds(org.Rn.eg.db, keys = genes, keytype = "SYMBOL",
                              column = "ENTREZID", multiVals = "first")
  if (mean(is.na(m1)) > 0.6) {
    m2 <- AnnotationDbi::mapIds(org.Rn.eg.db, keys = genes, keytype = "ALIAS",
                                column = "ENTREZID", multiVals = "first")
    idx <- which(is.na(m1) & !is.na(m2))
    m1[idx] <- m2[idx]
  }
  unname(stats::na.omit(m1))
}


amp <- apply(centroids, 2, function(v) diff(range(v, na.rm = TRUE)))
highlight <- names(sort(amp, decreasing = TRUE))[1:min(3, length(amp))]
highlight

for (clabel in highlight) {
  k <- as.integer(sub("^C","", clabel))
  genes_k <- dplyr::filter(cluster_df, Cluster == k) %>% dplyr::pull(Gene)   # SYMBOLs
  entrez  <- symbol_to_entrez(genes_k)
  message(sprintf("Cluster %d: %d genes, %d mapped to ENTREZ", k, length(genes_k), length(entrez)))
  if (length(entrez) < 10) next
  
  ego <- enrichGO(entrez, OrgDb = org.Rn.eg.db, keyType = "ENTREZID",
                  ont="BP", pAdjustMethod="BH", qvalueCutoff=0.05, readable=TRUE)
  
  ekegg <- enrichKEGG(entrez, organism = "rno",
                      pAdjustMethod="BH", qvalueCutoff=0.10)
  
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    p_go <- dotplot(ego, showCategory=10) + ggplot2::ggtitle(sprintf("GO BP - Cluster %d", k))
    ggsave(sprintf("final/GO_BP_cluster_%d.png", k), p_go, width=8, height=6, dpi=300)
    readr::write_csv(as.data.frame(ego), sprintf("final/GO_BP_cluster_%d.csv", k))
  }
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    p_kegg <- dotplot(ekegg, showCategory=10) + ggplot2::ggtitle(sprintf("KEGG - Cluster %d", k))
    ggsave(sprintf("final/KEGG_cluster_%d.png", k), p_kegg, width=8, height=6, dpi=300)
    # readable KEGG table (optional)
    try({
      ekegg_r <- setReadable(ekegg, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
      readr::write_csv(as.data.frame(ekegg_r), sprintf("final/KEGG_cluster_%d.csv", k))
    }, silent = TRUE)
  }
}

time_lookup <- pd_tissue %>%
  distinct(time_hours, time) %>%
  arrange(time_hours)

cent_df <- as.data.frame(centroids) %>%
  mutate(Sample = colnames(expr_z)) %>%
  left_join(pd_tissue %>% dplyr::select(Sample, time, time_hours), by = "Sample") %>%
  tidyr::pivot_longer(-c(Sample, time, time_hours), names_to = "Cluster", values_to = "Z") %>%
  group_by(Cluster, time, time_hours) %>%
  summarize(Z = mean(Z, na.rm = TRUE), .groups = "drop") %>%
  arrange(time_hours)


ggplot2::ggplot(cent_df, ggplot2::aes(time_hours, Z, group = Cluster)) +
  ggplot2::geom_line() + ggplot2::geom_point() +
  ggplot2::facet_wrap(~ Cluster, scales = "free_y") +
  ggplot2::scale_x_continuous(
    breaks = time_lookup$time_hours,
    labels = time_lookup$time
  ) +
  ggplot2::labs(
    title = paste("Cluster trajectories (", tissue_of_interest, ")", sep = ""),
    x = "Time", y = "Mean z-score"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

ggplot2::ggsave("final/kmeans_cluster_trajectories.png", width = 10, height = 7, dpi = 300)
