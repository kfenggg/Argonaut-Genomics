# TE contribution to HVG
library(Seurat)
library(dplyr)
library(ggplot2)

# read in preprocessed rds
combined.sobj <- readRDS("frontier2021_preprocessed.rds")

# read in te gtf and get unique
te.gtf <- rtracklayer::import('GRCh38_rmsk_TE.gtf')
te.gtf.df <- as.data.frame(te.gtf)
te.gene.ids <- te.gtf.df %>% select(gene_id) %>% unique()

# calc hvg and 
combined.sobj <- FindVariableFeatures(combined.sobj, 
                                      selection.method = "vst", 
                                      nfeatures = 10000)

top10000 <- head(VariableFeatures(combined.sobj), 10000)
top5000 <- head(VariableFeatures(combined.sobj), 5000)
top2000 <- head(VariableFeatures(combined.sobj), 2000)
top1000 <- head(VariableFeatures(combined.sobj), 1000)

sum(top10000 %in% te.gene.ids$gene_id)
sum(top5000 %in% te.gene.ids$gene_id)
sum(top2000 %in% te.gene.ids$gene_id)
sum(top1000 %in% te.gene.ids$gene_id)

# look at what they are
top1000[top1000 %in% te.gene.ids$gene_id]

top2000[top2000 %in% te.gene.ids$gene_id]
