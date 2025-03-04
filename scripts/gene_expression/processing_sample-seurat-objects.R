# Sample processing in Seurat

# utils
library(tidyverse)
library(Seurat)

# create seurat object from the gene count matrix
seurat <- CreateSeuratObject(counts = "matrix.tsv", min.cells = 3, min.features = 200)

# cell-level filtering
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
filtered_object <- subset(my_object, subset = 
                            nFeature_RNA > 0 & # min.feat
                            nFeature_RNA < 8000 & # max.feat
                            nCount_RNA > 0 & # min.count
                            nCount_RNA < 60000 & # max.count
                            percent.mt < 5 # mt
                          ) 

# normalization, vst, scaling, and pca 
seurat <- seurat %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

# identify PCs for clustering
ElbowPlot(seurat, ndims = 50)

# clustering
seurat <- seurat %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = seq(0.4, 2, by = 0.2)) %>%
  RunUMAP(dims = 1:30)

# UMAP visualization
DimPlot(seurat) 