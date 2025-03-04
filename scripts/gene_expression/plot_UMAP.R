# UMAP plotting function, updating on DimPlot

# utils
library(tidyverse)
library(Seurat)
library(paletteer)

# function
Plot_UMAP <- function(seurat_object, reduction = "umap", group = NULL, split = NULL, levels = NULL, colors = paletteer_d("ggthemes::Tableau_20"), 
                     pt_size = 0.01, label = FALSE, label_size = 8, text_size = 18, title = NULL, legend_title = NULL) {
  
  # order the plotting levels
  if (!is.null(levels)) {seurat_object@meta.data[[group]] <- factor(seurat_object@meta.data[[group]], levels = levels)}

  # plot
  DimPlot(object = seurat_object, reduction = reduction, group.by = group, split.by = split, cols = colors,# order = levels,
          pt.size = pt_size, label = label, label.size = label_size, repel = TRUE, raster = FALSE) +
    ggtitle(title) +
    xlab("") + 
    ylab("") +
    scale_color_manual(values = colors, name = legend_title) +
    theme_void() +
    theme(aspect.ratio = 1,
          text = element_text(size = text_size),
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = text_size),
          legend.position = "none",
          axis.title.y = element_text(hjust = 0.05, angle = 90), 
          axis.title.x = element_text(hjust = 0.05),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(1, 1, 1 ,1),
          panel.background = element_blank()
    )
}