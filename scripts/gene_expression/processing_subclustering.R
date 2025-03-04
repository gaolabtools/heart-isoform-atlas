# Subclustering processing and integration

# utils
library(tidyverse)
library(Seurat)

# subset each cell type
merged_seurat <- SetIdent(merged_seurat, value = "major_cell_type")
sub_seurat <- list("CM" = subset(merged_seurat, idents = "CM"),
                   "EC" = subset(merged_seurat, idents = "EC"),
                   "FB" = subset(merged_seurat, idents = "FB"),
                   "Mur" = subset(merged_seurat, idents = "Mur"),
                   "Mye" = subset(merged_seurat, idnets = "Mye"))

# process individual subset data
sub_seurat <- lapply(sub_seurat, function(x) {
  x %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA()
})

# integrate individual subset data
sub_seurat <- lapply(sub_seurat, function(x) {
  x %>%
    IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony") %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = 0.4, cluster.name = "harmony_clusters") %>%
    RunUMAP(reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")
})