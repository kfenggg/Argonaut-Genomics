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
featureCounts.outputs <- list.files(path = '2a/', pattern = '.rds', full.names = TRUE)
length(featureCounts.outputs)
# get file names instead of paths to make sure everything is in order
new.names <- list.files(path = '2a/', pattern = '.rds')
# drop .rds
new.names <- str_extract(new.names, pattern = "[^.rds]+")

all(cell.barcodes$old == new.names)


# load in our cellmetadata 
cell.metadata.complete <- read.csv('cell_metadata_complete.csv', row.names = NULL)
# drpo X column
cell.metadata.complete$X <- NULL
dim(cell.metadata.complete)
head(cell.metadata.complete, 5)

# number of barcodes we have that is in the cell metadata
table(cell.barcodes$new %in% cell.metadata.complete$barcode) # 1151 cells

# filter down to just the 1151 so we can run in 32gb of ram
cells.to.keep <- cell.barcodes$new %in% cell.metadata.complete$barcode
featureCounts.outputs <- featureCounts.outputs[cells.to.keep]
new.names <- cell.barcodes$new[cells.to.keep]


### CREATE COUNTS MATRIX 2a
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
saveRDS(counts.matrix, file = 'gse118389_2a_1151cells_counts_matrix.rds')

# filter our metadata down to just the cells we have and create a seurat object 
# subset the metadata to our 1151 cells
cell.metadata.subset <- cell.metadata.complete %>% filter(barcode %in% new.names)

# check the barcode order of counts matrix and cell metadata
table(cell.metadata.subset$barcode == colnames(counts.matrix), useNA = 'a')

# change none to Non-malignant and Normal to Normal-like
cell.metadata.subset$tumor.type[cell.metadata.subset$tumor.type == 'none'] <- 'Non-malignant'
cell.metadata.subset$tumor.type[cell.metadata.subset$tumor.type == 'Normal'] <- 'Normal-like'

head(cell.metadata.subset)
table(cell.metadata.subset$tumor.type)

rownames(cell.metadata.subset) <- cell.metadata.subset$barcode

# save cell metadata subset as rds
saveRDS(cell.metadata.subset, file = 'cell_metadata_subset_1151cells.rds')

##################################
### LOAD TO MAKE SEURAT OBJECT ###
##################################
# load in rds objects to make seurat object to save run time 
# min.features = 0 to keep all cells 
counts.matrix <- readRDS('gse118389_2a_1151cells_counts_matrix.rds')
cell.metadata.subset <- readRDS('cell_metadata_subset_1151cells.rds')
sobj.2a <- CreateSeuratObject(counts.matrix, 
                           project = 'GSE118389_2a', 
                           assay = 'TE_2a',
                           min.cells = 300, 
                           min.features = 0,
                           meta.data = cell.metadata.subset)
# EC2 saving
saveRDS(sobj.2a, file = 'gse118389_sobj2a_1151cells_300minCells.rds')
# local saving
saveRDS(sobj.2a, file = 'sobj_1151cells_rds/gse1181389_sobj2a_1151cells.rds')


### Check different sobj's
### LOAD sobj from EC2
sobj.500 <- readRDS('gse118389_sobj2a_1151cells_500minCells.rds')
sobj.300 <- readRDS('gse118389_sobj2a_1151cells_300minCells.rds')
sobj.300
