# Plot DEG volcano plot 

# utils
library(tidyverse)

# function
Plot_DEG_volcano <- function(degs, ct, nudge = 10, fc_x_limit = 8, text_size = 8) {
  
  p1 <- degs %>%
    # which cell type
    dplyr::filter(cell_type == ct) %>%
    # y p value limit for plotting
    mutate(p_val_adj = ifelse(p_val_adj <= 1e-250, 1e-250, p_val_adj)) %>%
    # ggplot
    ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    # points
    geom_point(color = "#737373", size = 1.6, pch = 16) +
    geom_point(data = . %>% dplyr::filter(avg_log2FC >= 0.6 & p_val_adj <= 0.05), color = "#F39B7F", size = 1, pch = 16) +
    geom_point(data = . %>% dplyr::filter(avg_log2FC <= -0.6 & p_val_adj <= 0.05), color = "#91D1C2", size = 1, pch = 16) +
    #
    ggtitle(ct) +
    coord_cartesian(clip = "off") +
    xlab(expression("log"[2] * "(fold change)")) +
    ylab(expression("-log"[10] * "(adj p-value)")) +
    scale_y_continuous(limits = c(0, 250), 
                       breaks = c(0, 50, 100, 150, 200, 250), 
                       labels = c("0", "50", "100", "150", "200", "\u2265 250")) +
    scale_x_continuous(limits = c(-fc_x_limit, fc_x_limit), 
                       breaks = c(-5, 0, 5)) +
    theme_bw() +
    theme(text = element_text(size = text_size),
          axis.text = element_text(size = text_size, color = "black"),
          plot.title = element_text(size = text_size, hjust = 0.5),
          aspect.ratio = 1)
  
  p1
}
