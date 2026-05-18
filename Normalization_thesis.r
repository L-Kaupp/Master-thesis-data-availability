#Install all packages and libraries

install.packages("SeuratObject")
install.packages("BiocManager")
BiocManager::install(c(
  "rhdf5",
  "HDF5Array",
  "zellkonverter"
))
BiocManager::install("glmGamPoi")
library(Seurat)
library(SingleCellExperiment)
library(SeuratObject)
library(zellkonverter)
#Define Working directory
setwd('/your/path')


#Load in AnnData object from QC script
sce <- readH5AD('/your/path/to/file/QC_AnnDataObject.h5ad', reader = "R")
sce 

# Convert SCE object to Seurat object
assayNames(sce)
counts_mat <- assay(sce, "X")
sobj <- CreateSeuratObject(
  counts = counts_mat
)

#Normalization (Log)
sobj <- NormalizeData(
  sobj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

#Normalization and variance stabilization (SCT)
options(future.globals.maxSize = 1000 * 1024^2)  
sobj <- SCTransform(sobj, variable.features.n = 1000)
seed.use = 1 

#Run PCS, UMAP, TSNE
sobj <- RunPCA(sobj, npcs = 50, verbose = FALSE)
ElbowPlot(sobj, reduction = "pca", ndims = 50)
sobj <- RunTSNE(sobj, dims = 1:12)
sobj <- RunUMAP(sobj, dims = 1:12)


#Clustering
# SNN Graph Construction
sobj <- FindNeighbors(sobj, reduction = "pca", k.param = 20, dims = 1:30)
# Cluster Determination:
sobj  <- FindClusters(sobj, resolution = seq(0.1, 0.9, by = 0.1), n.start = 100, n.iter = 100)

saveRDS(sobj, '/your/path') 
