library(tidyverse)
library(iso_seurat)
library(patchwork)

Plot_GeneandIsoformData <- function(iso_seurat, gene_seurat, gtf, gene, cell_type, group, groups, vln_fill, bar_colors, 
                             text_size = 8, UMI_threshold = 10) {
  
  # transcripts corresponding to genes
  tx <- gtf %>%
    dplyr::filter(gene_name %in% gene) %>%
    pull(transcript_ID)
  tx <- tx[tx %in% rownames(iso_seurat)]
  
  # create metadata column for aggregation
  iso_seurat$cluster_ID <- iso_seurat[[group]]
  
  # aggregate counts
  ps_data <- AggregateExpression(iso_seurat, group.by = "cluster_ID", features = tx, assays = "RNA")$RNA
  
  # prepare plotting data
  plot_data <- ps_data %>%
    as.data.frame() %>%
    # add gene IDs
    rownames_to_column(var = "transcript_ID") %>%
    left_join(gtf[c("transcript_ID", "gene_name", "transcript_name")], by = "transcript_ID") %>%
    # calculate isoform proportions
    pivot_longer(-c(transcript_name, transcript_ID, gene_name), names_to = "group", values_to = "count") %>%
    group_by(group) %>%
    mutate(pct = count / sum(count)) %>%
    ungroup() %>%
    # small pct as 'Other'
    mutate(transcript_name = ifelse(pct > 0.05, transcript_name, "Other")) %>%
    group_by(transcript_name, group) %>%
    mutate(pct = ifelse(transcript_name == "Other", sum(pct), pct)) %>%
    # plot levels, put Other as last
    mutate(transcript_name = factor(transcript_name),
           transcript_name = fct_relevel(transcript_name, "Other", after = Inf)) %>%
    # filter groups
    dplyr::filter(group %in% groups)
  
  # isoform ratios stacked bar  
  plot_data <- plot_data %>%
    distinct(transcript_name, pct) %>%
    ungroup() 
  
  # stacked bar plot
  p1 <- plot_data %>%
    # ggplot
    ggplot(aes(x = group, y = pct, fill = transcript_name)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
    #
    xlab("") +
    ylab("\nIsoform ratio") +
    scale_y_continuous(expand = c(0,0), 
                       limits = c(0, 1), 
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    coord_cartesian(clip = "off") +
    scale_fill_manual(values = bar_colors, name = NULL) +
    theme_bw() +
    theme(text = element_text(size = text_size),
          axis.text = element_text(size = text_size, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.text = element_text(size = text_size),
          legend.title = element_text(face = "italic"),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line.y.left = element_line(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(0,0,0,0)) 
  
  # gene expression violin plot
  # if the gene is found in the expression data
  if (gene %in% rownames(gene_seurat)) {
    # single cell expr values of the gene
    expr_data <- FetchData(gene_seurat, vars = gene) %>% 
      dplyr::rename("my_gene" = gene)
    # add cluster information
    expr_data <- merge(expr_data, gene_seurat[[group]], by = 0)
    colnames(expr_data)[3] <- "group"
    max_y <- ceiling(max(expr_data$my_gene) * 2) / 2 # max expression for plot limits to nearest 0.5
    
    p2 <- expr_data %>%
      mutate(group = str_replace(group, pattern = "_", "-")) %>%
      dplyr::filter(group %in% groups) %>%
      mutate(group = factor(group, levels = groups)) %>%
      # ggplot
      ggplot(aes(y = my_gene, x = group)) +
      geom_violin(linewidth = 0.2, fill = vln_fill, color = "black") +
      #
      scale_y_continuous(expand = c(0, 0), limits = c(0, max_y), breaks = c(seq(0, max_y, 1), max_y), labels = c(seq(0, max_y, 1), "")) +
      xlab(" ") + 
      ylab("Expression level") + 
      ggtitle(gene) +
      theme_bw() +
      theme(
        text = element_text(size = text_size),
        axis.text = element_text(color = "black", size = text_size),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(size = text_size, hjust = 0.5, face = "italic"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        plot.margin = margin(0,0,0)
      )
    
    patchwork::wrap_plots(p2, p1, ncol = 2, widths = c(0.8, 1), guides = "collect", axes = "collect_x")
    
  } else {
    # gene not found in expression data
    p1
  }
  
}
