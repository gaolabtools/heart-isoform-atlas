# Perform gene ontology overrepresentation analysis 

# utils
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)


# function for GOBP overrepresentation
run_GOBP_ORA <- function(degs) {
  res_list <- list()
  
  # perform the following for each cell type
  for (ct in unique(degs$cell_type)) {
    # significant positive DEGs only
    go.deg <- degs %>% 
      dplyr::filter(cell_type == ct & p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.6)
    
    # DEGs for non-diseased and heart-failure
    nd_genes <- bitr(go.deg %>% dplyr::filter(avg_log2FC < 0) %>% pull(gene), 
                     fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    hf_genes <- bitr(go.deg %>% dplyr::filter(avg_log2FC > 0) %>% pull(gene), 
                     fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    
    nd_res <- enrichGO(nd_genes$ENTREZID, 'org.Hs.eg.db', ont = "BP", minGSSize = 10, 
                       pvalueCutoff = 0.05, pAdjustMethod = "fdr", qvalueCutoff = 0.2, readable = TRUE)
    hf_res <- enrichGO(hf_genes$ENTREZID, 'org.Hs.eg.db', ont = "BP", minGSSize = 10, 
                       pvalueCutoff = 0.05, pAdjustMethod = "fdr", qvalueCutoff = 0.2, readable = TRUE)
    
    res_list[[ct]][["nd"]] <- nd_res
    res_list[[ct]][["hf"]] <- hf_res
    
  }
  res_list
}

# run on DEGs (HF vs ND)
results <- run_GOBP_ORA(x)
