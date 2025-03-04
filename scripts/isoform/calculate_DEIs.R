# Calculate differentially expressed isoforms

# utils
library(tidyverse)
library(Seurat)

# use FindAllMarkers() on isoform seurat
cm.dei <- FindAllMarkers(JoinLayers(subset(seurat, subset = subtype %in% c("CM1", "CM2", "CM3"))), assay = "RNA", only.pos = T)
ec.dei <- FindAllMarkers(JoinLayers(subset(seurat, subset = subtype %in% c("EC_cap", "EC_ven", "EC_art", "EC_cap-art", "EC_cap-ven", "EC_cap-ang"))), assay = "RNA", only.pos = T)
fb.dei <- FindAllMarkers(JoinLayers(subset(seurat, subset = subtype %in% c("FB1", "FB2", "FB3", "FB4"))), assay = "RNA", only.pos = T)
mur.dei <- FindAllMarkers(JoinLayers(subset(seurat, subset = subtype %in% c("PC", "SMC"))), assay = "RNA", only.pos = T)
mye.dei <- FindAllMarkers(JoinLayers(subset(seurat, subset = subtype %in% c("MP", "MP_Mo", "MP_prolif", "CD1C-DC", "CD1C+DC", "CD16+Mo", "Mast"))), assay = "RNA", only.pos = T)