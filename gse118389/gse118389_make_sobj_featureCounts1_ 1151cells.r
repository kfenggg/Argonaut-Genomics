# make sobj for featureCounts 1 
library(dplyr)
library(stringr)
library(Seurat)
library(RColorBrewer)
library(scuttle)
library(ggplot2)

#######################################
### CAN SKIP TO LOAD OBJECT SECTION ###
#######################################
# load cell barcodes that we fixed with Haiyin
cell.barcodes <- readRDS('df_fixed_cell_names.rds')
cell.barcodes
dim(cell.barcodes)

# read in each rds and construct a matrix
featureCounts.outputs <- list.files(path = '1/', pattern = '.rds', full.names = TRUE)
length(featureCounts.outputs)
# get file names instead of paths to make sure everything is in order
new.names <- list.files(path = '1/', pattern = '.rds')
# drop .rds
new.names <- str_extract(new.names, pattern = "[^.rds]+")

all(cell.barcodes$old == new.names)

# load subsetted metadata from making 
cell.metadata.complete <- readRDS('cell_metadata_subset_1151cells.rds')

# number of barcodes we have that is in the cell metadata
table(cell.barcodes$new %in% cell.metadata.complete$barcode) # 1151 cells

# filter down to just the 1151 so we can run in 32gb of ram
cells.to.keep <- cell.barcodes$new %in% cell.metadata.complete$barcode
featureCounts.outputs <- featureCounts.outputs[cells.to.keep]
new.names <- cell.barcodes$new[cells.to.keep]


### CREATE COUNTS MATRIX 1
# load first file and fix name
first.file <- readRDS(featureCounts.outputs[1])
colnames(first.file$counts) <- new.names[1]
head(first.file$counts)

# combine all the counts for all 1151 cells
counts.matrix <- as.data.frame(first.file$counts)
file.counter = 2 # to keep track of which file we are to rename columns
for (rds.path in featureCounts.outputs[2:length(featureCounts.outputs)]){
  print(file.counter)
  this.rds <- readRDS(rds.path)
  colnames(this.rds$counts) <- new.names[file.counter]
  counts.matrix[new.names[file.counter]] <- this.rds$counts
  file.counter = file.counter + 1 
}
# save the counts matrix as rds
saveRDS(counts.matrix, file = 'gse118389_1_1151cells_counts_matrix.rds')

# change none to Non-malignant and Normal to Normal-like
cell.metadata.complete$tumor.type[cell.metadata.complete$tumor.type == 'none'] <- 'Non-malignant'
cell.metadata.complete$tumor.type[cell.metadata.complete$tumor.type == 'Normal'] <- 'Normal-like'

# below is also used to fix counts matrix rownames
# since rownames are fixed we no longer need to use this for fixing DEA names
library(biomaRt)
# function to add a new column containing gene symbols to our FindMarkers df
add.gene.symbols <- function(marker.df){
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  genes <- rownames(marker.df)
  G_list <- getBM(filters = "ensembl_gene_id", 
                  attributes = c("ensembl_gene_id","hgnc_symbol"),
                  values = genes, mart = mart)
  print(all(rownames(marker.df) %in% G_list$ensembl_gene_id))
  gene.symbols <- c()
  for (id in rownames(marker.df)){
    gene.symbols <- c(gene.symbols, G_list[G_list$ensembl_gene_id == id,]$hgnc_symbol[1])
  }
  marker.df$gene.symbol <- gene.symbols
  return(marker.df)
}
# change ensemble id to gene id
test.matrix <- as.data.frame(rownames(counts.matrix))
colnames(test.matrix) <- 'ensembl.id'
rownames(test.matrix) <- test.matrix$ensembl.id

gene.id.matrix <- add.gene.symbols(test.matrix)
# fix empty spots
gene.id.matrix[gene.id.matrix$gene.symbol == '',]$gene.symbol <- gene.id.matrix[gene.id.matrix$gene.symbol == '',]$ensembl.id
rownames(counts.matrix) <- make.unique(gene.id.matrix$gene.symbol) # make unique adds '.1', '.2' etc to dups

# make sobj 
all(rownames(cell.metadata.complete) == colnames(counts.matrix))
sobj.1 <- CreateSeuratObject(counts.matrix, 
                              project = 'GSE118389_1', 
                              assay = 'RNA',
                              min.cells = 100, 
                              min.features = 1,
                              meta.data = cell.metadata.complete)
# save sobj
saveRDS(sobj.1, file = 'gse118389_featureCounts1_sobj_1151cells.rds')

### PREPROCESS 
##### SCTRANSFORM
# load 
sobj.1 <- readRDS('gse118389_featureCounts1_sobj_1151cells.rds')
# store mitochondrial percentage in object meta data
sobj.1 <- PercentageFeatureSet(sobj.1, pattern = "^MT-", col.name = "percent.mt")
summary(sobj.1$percent.mt)

# run sctransform
sobj.1 <- SCTransform(sobj.1, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
# These are now standard steps in the Seurat workflow for visualization and clustering
sobj.1 <- RunPCA(sobj.1, verbose = FALSE)
ElbowPlot(sobj.1)
sobj.1 <- RunUMAP(sobj.1, dims = 1:10, verbose = FALSE)

sobj.1 <- FindNeighbors(sobj.1, dims = 1:10, verbose = FALSE)
sobj.1 <- FindClusters(sobj.1, verbose = FALSE)
DimPlot(sobj.1, group.by = 'cell.type', cols = 'Set1')
DimPlot(sobj.1, group.by = 'tumor.type', cols = 'Set1')
#### NORMAL

#load
sobj.1 <- readRDS('gse118389_featureCounts1_sobj_1151cells.rds')
sobj.1 <- NormalizeData(sobj.1)
sobj.1 <- FindVariableFeatures(sobj.1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sobj.1)
sobj.1 <- ScaleData(sobj.1, features = all.genes)
sobj.1 <- RunPCA(sobj.1, features = VariableFeatures(object = sobj.1))
ElbowPlot(sobj.1, ndims = 30)

sobj.1 <- FindNeighbors(sobj.1, dims = 1:10)
sobj.1 <- FindClusters(sobj.1, resolution = 0.5)
sobj.1 <- RunUMAP(sobj.1, dims = 1:10)
DimPlot(sobj.1, reduction = "umap", group.by = 'cell.type', cols = 'Set1')
DimPlot(sobj.1, reduction = "umap", group.by = 'tumor.type', cols = 'Set1')

# Gene markers
# malig vs non malig epithelial gene markesr 
ge.malig.markers <- FindMarkers(sobj.1, only.pos = T,
                                group.by = 'tumor.type',
                                ident.1 = c('Basal', 'Her2', 'LumA', 'LumB', 'Normal-like'),
                                ident.2 = 'Non-malignant')
# sort by log fold change
ge.malig.markers.sorted <- ge.malig.markers[order(-ge.malig.markers$avg_log2FC), ]
# epithelial vs non-epithelial markers
ge.epi.markers <- FindMarkers(sobj.1, only.pos = T,
                              group.by = 'cell.type',
                              ident.1 = 'epithelial',
                              ident.2 = c('Bcell', 'endothelial', 'macrophage', 'stroma', 'Tcell', 'unk'))
# sort by log fold change
ge.epi.markers.sorted <- ge.epi.markers[order(-ge.epi.markers$avg_log2FC), ]


head(ge.malig.markers.sorted)
# convert ensemble id's to gene symbols
# ge.malig.markers.sorted.symbols <- add.gene.symbols(ge.malig.markers.sorted)
# ge.epi.markers.sorted.symbols <- add.gene.symbols(ge.epi.markers.sorted)
# head(ge.malig.markers.sorted.symbols)
# head(ge.epi.markers.sorted.symbols)

# save 
# write.csv(ge.malig.markers.sorted.symbols, file = 'featureCounts1_malig_markers_1151cells.csv')
# write.csv(ge.epi.markers.sorted.symbols, file = 'featureCounts1_epi_markers_1151cells.csv')
