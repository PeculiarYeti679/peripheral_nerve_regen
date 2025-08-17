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
#STUCK HERE TRYING TO GET ENTREZ ID
#######
#GO/KEGG Enrichment 
library(clusterProfiler)
library(org.Rn.eg.db)
library(GEOquery)
library(dplyr)
library(stringr)
library(tidyr)

#convert the symbols in the table to entrez id
symbol_to_entrez <- function(genes) {
  ids <- AnnotationDbi::mapIds(org.Rn.eg.db, keys=genes, keytype="SYMBOL",
                               column="ENTREZID", multiVals="first")
  unname(stats::na.omit(ids))
}

gse <- getGEO("GSE30165", GSEMatrix = FALSE) 

gpl_id <- names(GPLList(gse))
head(gpl_id)
#download GPL tagse#download GPL table
gpl <- getGEO(gpl_id)
#probe annotations
p_ann <- Table(gpl)
head(ann)
#get possible columns that have the symbols
possible_cols <- c("ID", "GENE", "GENE_SYMBOL")
sym_col <- possible_cols[possible_cols %in% colnames(ann)][1]
stopifnot(!is.null(sym_col))



#get the clisters with the highest amplitudes (top 2)
amp <- apply(centroids, 2, function(v) diff(range(v, na.rm=TRUE)))
highlight <- names(sort(amp, decreasing=TRUE))[1:min(3,length(amp))]

for (clabel in highlight) {

  k <- as.integer(sub("^C","", clabel))
  genes_k <- dplyr::filter(cluster_df, Cluster==k) |> dplyr::pull(Gene)
  entrez <- symbol_to_entrez(genes_k)
  if (length(entrez) < 10) next
  
  ego <- enrichGO(entrez, OrgDb=org.Rn.eg.db, keyType="ENTREZID",
                  ont="BP", pAdjustMethod="BH", qvalueCutoff=0.05, readable=TRUE)
  ekegg <- enrichKEGG(entrez, organism="rno", pAdjustMethod="BH", qvalueCutoff=0.10)
  
  if (!is.null(ego) && nrow(as.data.frame(ego))>0) {
    p <- dotplot(ego, showCategory=10) + ggtitle(paste("GO BP - Cluster", k))
    ggsave(paste0("final/GO_BP_cluster_", k, ".png"), p, width=8, height=6, dpi=300)
    readr::write_csv(as.data.frame(ego), paste0("final/GO_BP_cluster_", k, ".csv"))
  }
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg))>0) {
    p2 <- dotplot(ekegg, showCategory=10) + ggtitle(paste("KEGG - Cluster", k))
    ggsave(paste0("final/KEGG_cluster_", k, ".png"), p2, width=8, height=6, dpi=300)
    readr::write_csv(as.data.frame(ekegg), paste0("final/KEGG_cluster_", k, ".csv"))
  }
}
