# Single cell clustering using isoform count matrix

# utils
library(tidyverse)
library(Seurat)

# create seurat object from isoform count matrix
## isoform must be detected in at least 3 cells
seurat <- CreateSeuratObject(counts = count_matrix, names.delim = "_", meta.data = meta_data,
                              min.cells = 3, min.features = 0)

# normalization, vst, scaling, PCA
seurat <- seurat %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

# integration using the harmony method
seurat <- IntegrateLayers(object = seurat, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = 'harmony')
# clustering
seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:22)
seurat <- FindClusters(seurat, resolution = 0.4, cluster.name = "harmony_clusters")
seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:22, reduction.name = "umap.harmony")
