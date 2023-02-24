library(dplyr)
library(stringr)
library(Seurat)
library(RColorBrewer)
library(scuttle)
library(ggplot2)
library(ggrepel)

sobj.2a <- readRDS('gse118389_sobj2a_1151cells_300minCells.rds')
sobj.1 <- readRDS('gse118389_featureCounts1_sobj_1151cells.rds')
sobj.2a
sobj.1
#PREPROCESS TE sobj
sobj.2a <- NormalizeData(sobj.2a)
sobj.2a <- FindVariableFeatures(sobj.2a, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sobj.2a)
sobj.2a <- ScaleData(sobj.2a, features = all.genes)

# DEA
# TE markers
# malig vs non malig epithelial gene markers 
te.malig.markers <- FindMarkers(sobj.2a, only.pos = T,
                                group.by = 'tumor.type',
                                ident.1 = c('Basal', 'Her2', 'LumA', 'LumB', 'Normal-like'),
                                ident.2 = 'Non-malignant')
# sort by log fold change
te.malig.markers.sorted <- te.malig.markers[order(-te.malig.markers$avg_log2FC), ]
# epithelial vs non-epithelial markers
te.epi.markers <- FindMarkers(sobj.2a, only.pos = T,
                              group.by = 'cell.type',
                              ident.1 = 'epithelial',
                              ident.2 = c('Bcell', 'endothelial', 'macrophage', 'stroma', 'Tcell', 'unk'))
# sort by log fold chante
te.epi.markers.sorted <- te.epi.markers[order(-te.epi.markers$avg_log2FC), ]


head(te.malig.markers.sorted)
head(te.epi.markers.sorted)

# SAVE
write.csv(te.malig.markers.sorted, file = 'featureCounts2a_malig_markers.csv')
write.csv(te.epi.markers.sorted, file = 'featureCounts2a_epithelial_markers.csv')

### COMBINE SOBJS
### ADD GENE EXPRESSION COUNTS MATRIX AS ASSAY TO SOBJ
all(Cells(sobj.1) == Cells(sobj.2a)) # check cell names are the same

# create new assay to store ge matrix
ge.assay <- GetAssayData(object = sobj.1, slot = "counts")
ge.assay <- CreateAssayObject(counts = ge.assay)
# add to sobj
sobj <- sobj.2a
sobj[["RNA"]] <- ge.assay
DefaultAssay(sobj) <- "RNA"

############################
### SCTRANSFORM ON GENES ###
############################
# store mitochondrial percentage in object meta data
sobj <- PercentageFeatureSet(sobj, pattern = "^MT-", col.name = "percent.mt")
summary(sobj$percent.mt)

# run sctransform
sobj <- SCTransform(sobj, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
# These are now standard steps in the Seurat workflow for visualization and clustering
sobj <- RunPCA(sobj, verbose = FALSE)
ElbowPlot(sobj)
sobj <- RunUMAP(sobj, dims = 1:30, verbose = FALSE)

sobj <- FindNeighbors(sobj, dims = 1:30, verbose = FALSE)
sobj <- FindClusters(sobj, verbose = FALSE)

DimPlot(sobj, group.by = 'cell.type', cols = 'Set1')
DimPlot(sobj, group.by = 'tumor.type', cols = c('Pink', 'Green', 'Blue', 'Purple', 'Black', 'Red'))

# add malig cell.type classifier to the meta data
malig.status <- c()
for (idx in 1:length(Cells(sobj))) {
  if((sobj@meta.data$cell.type[idx] == 'epithelial') & (sobj@meta.data$tumor.type[idx] == 'Non-malignant')){
    malig.status <- c(malig.status, 'Non-malignant-epithelial')
  }
  else if(sobj@meta.data$cell.type[idx] == 'epithelial'){
    malig.status <- c(malig.status, 'Malignant-epithelial')
  }
  else{
    malig.status <- c(malig.status, sobj@meta.data$cell.type[idx])
  }
}
sobj@meta.data$malig.cell.type <- malig.status
# replot umap with this classifier 
DimPlot(sobj, group.by = 'malig.cell.type', cols = c('Cyan', 'Blue', 'Green', 'Red', 'Black', 'Yellow', 'Purple', 'Pink'))

# feature plots of top 10 TE's 
FeaturePlot(sobj, reduction = 'umap', features = rownames(te.malig.markers.sorted)[4:6])
FeaturePlot(sobj, reduction = 'umap', features = rownames(te.epi.markers.sorted)[1:10])


# read in te gtf
te.gtf <- rtracklayer::import('../GRCh38_rmsk_TE.gtf')
te.gtf.df <- as.data.frame(te.gtf)
# for features you want to look up make sure you convert '-' to '_'
te.gtf.df[te.gtf.df$transcript_id == 'AluSx1_dup38012',]
te.gtf.df[te.gtf.df$transcript_id == 'L1MB7_dup1019',]

# get width's for the top 10 TE's
top.10.te.malig <- rownames(te.malig.markers.sorted)[1:10]
top.10.te.epi <- rownames(te.epi.markers.sorted)[1:10]

# function output TE + width 
get.width <- function(te.list){
  cat('name chr start end strand width', '\n')
  for (te in te.list){
    
    cat(te, te.gtf.df[te.gtf.df$transcript_id == gsub('-', '_', te),]$seqnames,
        te.gtf.df[te.gtf.df$transcript_id == gsub('-', '_', te),]$start,
        te.gtf.df[te.gtf.df$transcript_id == gsub('-', '_', te),]$end,
        te.gtf.df[te.gtf.df$transcript_id == gsub('-', '_', te),]$strand,
        te.gtf.df[te.gtf.df$transcript_id == gsub('-', '_', te),]$width, '\n')
  }
}
get.width(top.10.te.malig)
get.width(top.10.te.epi)
################################
### DEA FOR GENES, SCT ASSAY ###
################################
# Gene markers
# malig vs non malig epithelial gene markers 
gene.malig.markers <- FindMarkers(sobj, only.pos = T,
                                group.by = 'tumor.type',
                                ident.1 = c('Basal', 'Her2', 'LumA', 'LumB', 'Normal-like'),
                                ident.2 = 'Non-malignant')
# sort by log fold change
gene.malig.markers.sorted <- gene.malig.markers[order(-gene.malig.markers$avg_log2FC), ]
# epithelial vs non-epithelial markers
gene.epi.markers <- FindMarkers(sobj, only.pos = T,
                              group.by = 'cell.type',
                              ident.1 = 'epithelial',
                              ident.2 = c('Bcell', 'endothelial', 'macrophage', 'stroma', 'Tcell', 'unk'))
# sort by log fold chante
gene.epi.markers.sorted <- gene.epi.markers[order(-gene.epi.markers$avg_log2FC), ]


head(gene.malig.markers.sorted)
head(gene.epi.markers.sorted)

# feature plots of top 10 TE's 
FeaturePlot(sobj, reduction = 'umap', features = rownames(gene.malig.markers.sorted)[1:10])
FeaturePlot(sobj, reduction = 'umap', features = rownames(gene.epi.markers.sorted)[1:10])

# SAVE 
write.csv(gene.malig.markers.sorted, file = 'featureCounts1_malig_markers_1151cells.csv')
write.csv(gene.epi.markers.sorted, file = 'featureCounts1_epi_markers_1151cells.csv')

# save the combined sobj
saveRDS(sobj, 'gse118389_combined_sobj_1_2a_1oct22.rds')


#### Load in 
sobj <- readRDS('gse118389_combined_sobj_1_2a_1oct22.rds')
# change unk to unknow in metadata
sobj@meta.data[sobj@meta.data$malig.cell.type == 'unk',]$malig.cell.type <- 'unknown'

# read in te csvs
te.malig.markers.sorted <- read.csv(file = 'featureCounts2a_malig_markers.csv', row.names = 1)
te.epi.markers.sorted <- read.csv(file = 'featureCounts2a_epithelial_markers.csv', row.names = 1)

te.gtf.df[te.gtf.df$transcript_id == 'Charlie6_dup220',]

## gene set score for L1MB8 dups in malig

DefaultAssay(sobj) <- 'TE_2a'
sobj <- AddModuleScore(sobj, features = list(top.10.te.malig[2:6]), name = 'L1MB_dup_score')



FeaturePlot(sobj, features = 'L1MB_dup_score1')

# notch expression, gene close to the L1MB's
DefaultAssay(sobj) <- 'SCT'
FeaturePlot(sobj, features = c('NOTCH2NLB', 'NOTCH2NLA', 'NOTCH2NLC',
                               'NBPF26', 'NOTCH2'))
