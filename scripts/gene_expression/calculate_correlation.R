# Calculate gene expression correlation between in-house and public data

# utils
library(tidyverse)
library(Seurat)

# calculate average expression
avg_data <- AverageExpression(large_seurat, group.by = c("source", "Condition", "cell_type"), 
                               layer = "data", assays = "RNA")

# prepare data
avg_data <- avg_data$RNA %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_ID") %>%
  pivot_longer(-gene_ID) %>%
  separate(name, into = c("source", "condition", "cell_type"), sep = "_") %>%
  pivot_wider(id_cols = c("gene_ID", "condition", "cell_type"), names_from = "source", values_from = "value")

# calculate Spearman coefficients for each cell type under each condition
avg_data %>%
  group_by(condition, cell_type) %>%
  mutate(corr_coef = cor(`Northwestern`, `Broad_MGH`, method = "spearman"),
         corr_coef = round(corr_coef, 3),
         corr_coef = if_else(row_number() == 1, corr_coef, NA_real_)) 