# perform GSEA on DEGs

# utils
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)


run_GSEA <- function(degs, ct = NULL) {
  # get and rank all degs and FC for either bulk or cell type
  if (!is.null(ct)) {
    ct_data <- degs %>%
      dplyr::filter(cell_type == ct) %>%
      dplyr::select(gene, avg_log2FC) %>%
      arrange(desc(avg_log2FC)) %>%
      filter(!is.na(avg_log2FC))
  } else {
    ct_data <- degs %>%
      rownames_to_column(var = "gene") %>%
      dplyr::select(gene, avg_log2FC) %>%
      arrange(desc(avg_log2FC)) %>%
      filter(!is.na(avg_log2FC))
  }
 
  # get entrez ids
  entrez <- bitr(ct_data %>% pull(gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") %>%
    mutate(ENTREZID = as.numeric(ENTREZID)) %>%
    dplyr::rename("gene" = SYMBOL) 
  
  # sorted named vector 
  gsea_data <- ct_data %>%
    right_join(., entrez, by = "gene") %>%
    dplyr::select(-gene) %>%
    arrange(desc(avg_log2FC)) %>%
    dplyr::select(ENTREZID, avg_log2FC) %>%
    deframe()
  
  # get entrez ids in gene sets
  t2g <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>%
    dplyr::select(gs_name, entrez_gene)

  # run gsea
  GSEA(gsea_data, TERM2GENE = t2g)
}

all.gsea <- runGSEA(bulkdeg)
cm.gsea <- runGSEA(deg, ct = "Cardiomyocyte")
ec.gsea <- runGSEA(deg, ct = "Endothelial")
fb.gsea <- runGSEA(deg, ct = "Fibroblast")
smc.gsea <- runGSEA(deg, ct = "Smooth muscle")
pc.gsea <- runGSEA(deg, ct = "Pericyte")
mye.gsea <- runGSEA(deg, ct = "Myeloid")
lec.gsea <- runGSEA(deg, ct = "Lymphoid")
lym.gsea <- runGSEA(deg, ct = "Lymphatic")
edc.gsea <- runGSEA(deg, ct = "Endocardial")
neu.gsea <- runGSEA(deg, ct = "Neuronal")


# visualization
gseaplot2(all.gsea, geneSetID = c(3, 20, 26, 16, 19, 45), subplots = 1:2, pvalue_table = FALSE) 