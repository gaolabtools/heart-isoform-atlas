# Calculate differentially expressed genes (DEGs) 

# utils
library(Seurat)

# calculate DEGs for each major cell type
merged_seurat <- SetIdent(merged_seurat, value = "major_cell_type")
DEG_celltype <- FindAllMarkers(JoinLayers(merged_seurat), only.pos = TRUE)

# calculate DEGs (HF vs ND) for each major cell type
DEG_HFvsND <- list()
for (celltype in unique(merged_seurat$cell_type_integrated)) {
  deg_list[[celltype]] <- 
    FindMarkers(JoinLayers(subset(merged_seurat, idents = celltype)), group.by = "Condition", ident.1 = "HF", ident.2 = "ND", only.pos = FALSE) %>%
    mutate(cell_type = celltype)
}
DEG_HFvsND <- DEG_HFvsND %>%
  purrr::reduce(., rbind)

# calculate DEGs (HF vs ND) bulk
bulkdeg <- FindMarkers(JoinLayers(merged_seurat), group.by = "Condition", ident.1 = "HF", ident.2 = "ND")

# calculate DEGs for each subcluster
DEG_subcluster <- lapply(sub_seurat, function(x) {
  FindAllMarkers(JoinLayers(x), only.pos = TRUE)
})
