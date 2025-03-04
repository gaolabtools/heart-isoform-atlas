# Integration of gene expression data using Seurat and Harmony

# utils
library(tidyverse)
library(Seurat)

# read in seurat objects from each sample
seurat_list <- lapply(Sys.glob("*_seurat.rds"), readRDS)

# merge seurat objects
merged_seurat <- lapply(seurat_list, function(x) {merge(x[[1]], y = x[-1], merge.data = FALSE)})

# clustering without integration
merged_seurat <- merged_seurat %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:50) %>%
  FindClusters(resolution = seq(0.4, 2, by = 0.2)) %>%
  RunUMAP(dims = 1:50) 

# harmony integration and clustering
merged_seurat <- IntegrateLayers(merged_seurat, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony")
merged_seurat <- merged_seurat %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.4, cluster.name = "harmony_clusters") %>%
  RunUMAP(reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
DimPlot(merged_seurat, group.by = c("harmony_clusters", "orig.ident"), reduciton = "umap.harmony", label = FALSE)


## integration of in-house and public data
# create metadata column
inhouse_seurat <- merged_seurat
inhouse_seurat[["source"]] <- "Northwestern"
public_seurat[["source"]] <- "Broad_MGH"

large_seurat <- merge(inhouse_seurat, public_seurat, merge.data = F)

# process 
large_seurat <- large_seurat %>%
  NormalizeData(assay = "RNA") %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

ElbowPlot(large_seurat) 

# clustering
large_seurat <- large_seurat %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.6) %>%
  RunUMAP(dims = 1:30)

# harmony integration
large_seurat <- Harmony_integration(large_seurat, sigma = 0.2)
large_seurat <- Harmony_clustering(large_seurat, dims = 30, res = 0.5)
