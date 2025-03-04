# Plot isoform tracks

# utils
library(tidyverse)
library(ggtranscript)
library(paletteer)

# function
Plot_IsoTrack <- function(gtf, gene, transcripts, text_size = 18) {
  
  # filter for gene and transcripts 
  annotation <- gtf %>% 
    dplyr::filter(gene_name == gene) %>%
    dplyr::select(seqnames, start, end, strand, type, gene_name, transcript_type, transcript_id, transcript_name)%>% 
      dplyr::filter(transcript_name %in% transcripts)
  
  # exons
  exons <- annotation %>% 
    dplyr::filter(type == "exon")
  
  # scale coordinates
  rescaled <- shorten_gaps(
    exons = exons, 
    introns = to_intron(exons, "transcript_id"), 
    group_var = "transcript_id")
  
  p1 <- rescaled %>%
    dplyr::filter(type == "exon") %>%
    # ggplot
    ggplot(aes(xstart = start, xend = end, y = transcript_name)) +
    geom_range(aes(fill = transcript_type), color = "transparent") +
    geom_intron(data = rescaled %>% dplyr::filter(type == "intron"), color = "grey40", arrow.min.intron.length = 200) +
    coord_cartesian(clip = "off") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = text_size),
          axis.text = element_text(size = text_size, color = "black"),
          legend.text = element_text(size = text_size),
          plot.title = element_text(size = text_size, face = "italic"),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          plot.margin = margin(0, 0, 0, 0))
  
}