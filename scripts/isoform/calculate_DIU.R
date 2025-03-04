# Perform the Chi-sq test using aggregate UMI isoform counts

# utils
library(tidyverse)
library(Seurat)

# function
RunX2 <- function(seurat, gtf, DR_data, group.1, group.2, UMI_thresh = 10, sig_dpct_thresh = 0.1, sig_pval_thresh = 0.05 ) {
  
  # aggregate data (pseudobulk)
  ps_data <- AggregateExpression(seurat, group.by = "isoform_cluster", assays = "RNA")$RNA
  
  # prepare data
  ps_data <- ps_data %>%
    as.data.frame() %>%
    # add gene names
    rownames_to_column(var = "transcript_ID") %>%
    left_join(gtf[c("transcript_ID", "transcript_name", "gene_name")], by = "transcript_ID") %>%
    # long format
    pivot_longer(-c("transcript_ID", "transcript_name", "gene_name"), values_to = "count", names_to = "cluster")
  
  # determine genes to test using detection rate
  genes <- DR_data %>%
    # group1 and group2 filter
    dplyr::filter(cluster %in% c(group.1, group.2)) %>%
    # DR filter (gene is detected in at least 5% cells in both groups)
    group_by(gene_name) %>%
    dplyr::filter(all(DR >= 0.05)) %>%
    pull(gene_name) %>%
    unique()
  
  # counts data for sample and groups
  counts_data <-
    ps_data %>% 
    # group1 and group2 filter
    dplyr::filter(cluster %in% c(group.1, group.2)) %>%
    # genes to test
    dplyr::filter(gene_name %in% genes) %>%
    # wide format
    pivot_wider(id_cols = c(transcript_ID, transcript_name, gene_name), values_from = "count", names_from = cluster) %>%
    # UMI filter
    dplyr::filter(!if_all(where(is.numeric), ~ . <= UMI_thresh)) %>%
    # remove genes with only 1 tx
    group_by(gene_name) %>%
    dplyr::filter(n() > 1) %>%
    ungroup() %>%
    # rename group columns as counts
    dplyr::rename(cts.1 = group.1,
                  cts.2 = group.2) %>%
    # calculate difference in proportions
    group_by(gene_name) %>%  
    arrange(gene_name) %>%
    mutate(pct.1 = cts.1 / sum(cts.1),
           pct.2 = cts.2 / sum(cts.2),
           dpct = pct.1 - pct.2) %>%
    ungroup() %>%
    # add comparison column
    mutate(comparison = paste(group.1, group.2, sep = "vs"))
  
  # number of tests for pval correction
  n.tests <- length(unique(counts_data$gene_name))
  
  # run chi-square test
  chisq_res <- counts_data %>%
    group_by(gene_name) %>%
    mutate(chisq_pvalue = chisq.test(matrix(c(cts.1, cts.2), ncol = 2, byrow = FALSE))$p.value,
           chisq_stat = chisq.test(matrix(c(cts.1, cts.2), ncol = 2, byrow = FALSE))$statistic) %>%
    ungroup()
  
  # calculate adjusted pval
  chisq_res <- 
    chisq_res %>% 
    dplyr::select(gene_name, comparison, chisq_pvalue) %>%
    distinct() %>%
    mutate(adj_p_val = p.adjust(chisq_pvalue, method = "BH")) %>%
    left_join(chisq_res, by = c("gene_name", "comparison", "chisq_pvalue")) 
}