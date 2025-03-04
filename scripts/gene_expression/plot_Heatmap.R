# Plot heatmap of DEG expression

# utils
library(tidyverse)
library(Seurat)
library(ComplexHeatmap)
library(paletteer)

# function
Plot_Heatmap <- function (seurat_object, degs, gene_labels,
                                   celltype_column, condition_column, condition_colors, condition_legend_title,
                                   text_size = 8, heat_colors = NULL, colors_legend_title) {
  
  # calculate average expression from data slot, convert to matrix
  mydata <- AverageExpression(object = seurat_object, features = degs$gene, assays = "RNA", layer = "data", 
                              group.by = c(celltype_column, condition_column))
  mydata <- as.matrix(unlist(mydata$RNA))
  
  # reorder rows according to input gene list
  mydata <- mydata[degs$gene, ]

  # column
  celltype_column_vec <- str_split_i(colnames(mydata), pattern = "_", i = 1)
  
  # row colors and split by celltype
  celltype_row_vec <- degs$cluster
  
  # column colors of condition
  condition_column_vec <- str_split_i(colnames(mydata), pattern = "_", i = 2)
  condition_colors <- setNames(as.vector(condition_colors), condition_levels)

  # gene labels
  ha = rowAnnotation(foo = anno_mark(at = which(rownames(mydata) %in% gene_labels),
                                     labels = rownames(mydata)[which(rownames(mydata) %in% gene_labels)],
                                     labels_gp = gpar(fontsize = round(text_size * 0.8), fontface = "italic")))
  # plot
  Heatmap(
    # scaled matrix
    matrix = t(scale(t(mydata))),
    # splits
    row_split = celltype_row_vec, 
    column_split = celltype_column_vec,
    # customization
    if (!is.null(heat_colors)) {col = heat_colors}, 
    border = TRUE, 
    cluster_columns = FALSE, 
    cluster_rows = FALSE,
    show_column_names = FALSE, 
    show_row_names = FALSE,
    column_title_rot = 45, 
    column_title_gp = gpar(fontsize = text_size),
    row_title = NULL,
    # annotations
    right_annotation = ha,
    top_annotation = HeatmapAnnotation(cnd = condition_column_vec,
                                       col = list(cnd = condition_colors),
                                       show_annotation_name = c(F,F),
                                       annotation_label = condition_legend_title,
                                       simple_anno_size = unit(3, "mm"),
                                       annotation_legend_param = list(title_gp = gpar(fontsize = text_size, fontface = "bold"),
                                                                      labels_gp = gpar(fontsize = text_size))),
    heatmap_legend_param = list(title = colors_legend_title, direction = "vertical",
                                title_gp = gpar(fontsize = text_size, fontface = "bold"), 
                                labels_gp = gpar(fontsize = text_size),
                                legend_height = unit(3, "cm"), legend_width = unit(4, "cm"))
  )
}
