#Design matrices and analyze data

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

needed <- c(
  "AnnotationDbi", "org.Rn.eg.db", "TCseq", "clusterProfiler",
  "GEOquery", "limma", "janitor", "pheatmap", "ggplot2", "tidyverse", "ggrepel"
)

to_install <- setdiff(needed, rownames(installed.packages()))
if (length(to_install)) {
  BiocManager::install(to_install, ask = FALSE, update = FALSE)
}


update.packages(ask = FALSE, checkBuilt = TRUE)
install.packages("rmarkdown", dependencies = TRUE)

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Rn.eg.db)     # rat annotation
  library(TCseq)
  library(clusterProfiler)
  library(GEOquery)
  library(limma)
  
  library(janitor)
  library(pheatmap)
  library(ggplot2)
  library(tidyverse)
})
#read in saved data
#this is processed in 01_Data_Acquistion_Normalization
expr <- readRDS("data/processed/expression_matrix_processed.rds")
pheno <- read.csv("data/meta/phenodata.csv", stringsAsFactors = FALSE)

#check for name matching in columns
stopifnot(all(colnames(expr) == pheno$Sample))

#converting to variables to use as catergory predictors
pheno$tissue <- factor(pheno$tissue)  
pheno$time   <- factor(pheno$time)
#remove unused varaibles
pheno$tissue <- droplevels(pheno$tissue)
pheno$time   <- droplevels(pheno$time)

#make different groups for different times
design <- model.matrix(~ 0 + tissue:time, data = pheno)
colnames(design) <- make.names(colnames(design))

head(design)

#compare the data between 0d and 1d of type DRG
#or 0d and 1d of type SN
contrast_matrix <- makeContrasts(
  DRG_1d_vs_0d = tissueDRG.time1d - tissueDRG.time0d,
  SN_1d_vs_0d  = tissueSN.time1d  - tissueSN.time0d,
  levels = design
)
head(contrast_matrix)



#create a linear model
#fit linear model

fit <- lmFit(expr, design)

fit2 <- contrasts.fit(fit, contrast_matrix)
#ebayes used to help with variance in genes
fit2 <- eBayes(fit2)

topTable_DRG <- topTable(fit2, coef = "DRG_1d_vs_0d", number = Inf, adjust = "fdr")
head(topTable_DRG)


gpl <- GEOquery::getGEO("GPL7294", AnnotGPL = TRUE)
gpl_tbl <- Table(gpl) |> as_tibble() |> clean_names()
names(gpl_tbl)
#mape the probeID to the gene aymbol, entrez, or gene
annot <- gpl_tbl %>%
  transmute(
    ProbeID = id,                            
    SYMBOL  = na_if(gene_symbol, ""),         
    ENTREZ  = ifelse(is.na(gene), NA_character_, as.character(gene)), 
    GENE    = coalesce(gene_name, description) 
  ) %>%
  distinct(ProbeID, .keep_all = TRUE)
#add the annotation to the topTable
topTable_DRG_annot <- topTable_DRG %>%
  rownames_to_column("ProbeID") %>%
  left_join(annot, by = "ProbeID")

#create volcano plot
library(ggplot2)

#red = significant
#grey = not sig
#used to spot strong DEGs (Differential Expressed Genes)
ggplot(topTable_DRG_annot, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = adj.P.Val < 0.05 & abs(logFC) >= 1)) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: DRG 1d vs 0d", x = "log2 Fold Change", y = "-log10 adj P-value")


#heatmap
library(pheatmap)

#get the top 50 probes
probes <- intersect(topTable_DRG_annot$ProbeID, rownames(expr))
top50  <- head(probes, 50)
#check to make sure all data is there
stopifnot(all(colnames(expr) == pheno$Sample))
#build the matrix
heat_mat <- expr[top50, , drop = FALSE]

#build ann which will keep tissue, sample and time 
#no duplicates
#ann is cleaned metadata
ann <- pheno %>%
  transmute(Sample = Sample,
            tissue = ifelse(is.na(tissue), "Unknown", tissue),
            time   = ifelse(is.na(time),   "Unknown", time)) %>%
  distinct(Sample, .keep_all = TRUE) %>%
  column_to_rownames("Sample")

#convert to factors
#remove empty level
ann$tissue <- droplevels(factor(ann$tissue))
ann$time   <- droplevels(factor(ann$time))

#ensure columns match
stopifnot(identical(colnames(heat_mat), rownames(ann)))

pheatmap(
  heat_mat,
  scale = "row",
  annotation_col = ann,
  show_rownames = FALSE,
  show_colnames = FALSE
)

#get sig DEG
topTable_DRG_annot <- topTable_DRG_annot %>%
  mutate(sig = adj.P.Val < 0.05 & abs(logFC) >= 1)

#mark sig (red) 
#not sig grey 
ggplot(topTable_DRG_annot, aes(logFC, -log10(adj.P.Val), color = sig)) +
  geom_point(alpha = 0.7, size = 1.2) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red"),
                     labels = c("Not sig", "adj.P.Val < 0.05 & |log2FC| ≥ 1"),
                     name = "DEG") +
  theme_minimal() +
  labs(title = "DRG 1d vs 0d", x = "log2 Fold Change", y = "-log10 adj P-value")

#save data
write.csv(topTable_DRG_annot, "outputs/tables/DRG_1d_vs_0d_DEGs_annotated.csv", row.names = FALSE)
saveRDS(topTable_DRG_annot,  "outputs/tables/DRG_1d_vs_0d_DEGs_annotated.rds")
