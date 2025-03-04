# Relationship between transcript diversity and isoform proportion

# utils
library(tidyverse)
library(Seurat)

# aggregate counts
agg_data <- AggregateExpression(iso_sobject, group.by = "isoform_cluster")$RNA %>%
  as.data.frame() %>%
  rownames_to_column(var = "transcript_ID") %>%
  left_join(gtf[c("transcript_ID", "gene_name")], by = "transcript_ID") %>%
  pivot_longer(-c(gene_name, transcript_ID), names_to = "celltype", values_to = "count")


calc.shannon <- function(x) { 
  top2 <- head(sort(x, decreasing = TRUE), 2)
  -sum(top2[top2 > 0] * log(top2[top2 > 0]))
}

# calculate diversity
div_data <- agg_data %>%
  dplyr::filter(count != 0) %>%
  group_by(gene_name, celltype) %>%
  dplyr::filter(n() != 1) %>%
  mutate(pct = count / sum(count),
         H = calc.shannon(pct),
  ) %>%
  ungroup() %>%
  arrange(celltype, gene_name)

# get only the proportion of the most dominant isoform
plot_data <- div_data %>%
  group_by(gene_name, celltype) %>%
  slice_max(pct, with_ties = FALSE) %>%
  ungroup()

# ggplot
plot_data %>%
  ggplot(aes(x = pct, y = H)) +
  geom_point(pch = 16, size = 0.5, color = "grey70") +
  geom_smooth(method = "loess", color = "#990F0F", linewidth = 0.7) +
  geom_hline(yintercept = 0.5, color = "#DF8F44", linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = 0, y = 0.5, vjust = 1.3, hjust = 0, fontface = "italic", label = "monoform", color = "#DF8F44", size = 18/.pt) +
  annotate("text", x = 0, y = 0.5, vjust = -0.3, hjust = 0, fontface = "italic", label = "polyform", color = "#DF8F44", size = 18/.pt) +
  scale_y_continuous(breaks = seq(0, 0.7, 0.1)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  coord_cartesian(clip = "off") +
  xlab("Dominant isoform prop") +
  ylab(expression("Diversity (" * italic(D)[italic(I)] * ")")) +
  theme_bw() +
  theme(text = element_text(size = 18),
        axis.text = element_text(size = 18, color = "black"), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 1)
