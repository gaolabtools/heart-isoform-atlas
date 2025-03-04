# Code for the analysis of isoform data

# utils
library(tidyverse)

## transcript filtering
# load metadata from gene expression
metadata <- read_tsv("metadata.tsv", col_names = TRUE, show_col_types = FALSE)

# load SQANTI output
sqanti <- read_delim("OUT.extended_annotation_classification.txt", col_names = TRUE, show_col_types = FALSE)

# FSM and nonFSM isoforms
FSM_tx <- sqanti %>%
  dplyr::filter(structural_category == "full-splice_match") %>%
  pull(isoform)

nonFSM_tx <- sqanti %>%
  dplyr::filter(structural_category != "full-splice_match") %>%
  dplyr::filter(within_CAGE_peak == TRUE) %>% # with 5' CAGE support
  pull(isoform)

# load isoform count matrix
count_mat <- read_tsv("matrix_isoform_all_Samples.tsv", 
                      col_select = c("#feature_id", metadata$cell_barcode) , col_names = TRUE, show_col_types = FALSE)
colnames(count_mat)[1] <- "transcript_ID"

# filter count matrix
filt_count_mat <- count_mat %>%
  dplyr::filter(transcript_ID %in% c(FSM_tx, nonFSM_tx))

# calculate detection rates (DR)
samples <- unique(metadata$orig.ident)
DR_list <- list()
for (sample in samples) {
  # cell barcodes for each sample
  cell_bc <-
    metadata %>%
    dplyr::filter(orig.ident == sample) %>%
    pull(cell_barcode)
  # count mat of sample
  DR_list[[sample]] <-
    filt_count_mat %>%
    dplyr::select("transcript_ID", cell_bc) %>%
    # calculate total UMIs for each transcript
    dplyr::mutate("{sample}" := rowSums(across(where(is.numeric)))) %>%
    dplyr::select("transcript_ID", sample)
}

DR_data <- DR_list %>%
  purrr::reduce(full_join, by = "transcript_ID")

# add known vs novel identifiers
DR_data <-
  sqanti %>%
  rename(transcript_ID = "isoform") %>%
  dplyr::select(transcript_ID, structural_category) %>%
  right_join(DR_data, by = "transcript_ID")

# determine how many samples each transcript is detected in
DR_data <-
  DR_data %>%
    mutate(n_UMI = rowSums(across(starts_with("HR"))),
           n_samples = rowSums(across(starts_with("HR"), ~ . > 0)))

# transcript filtering by UMI
# for nonFSM - found in at least 2 samples
tx_to_keep <- DR_data %>%
  dplyr::filter(n_UMI > 0) %>% # from removed cells
  dplyr::filter(structural_category == "full-splice_match" | (structural_category != "full-splice_match" & n_samples >= 2 & n_UMI >= 3)) %>%
  pull(transcript_ID)

umi_filt_count_mat <- 
  filt_count_mat %>%
  dplyr::filter(transcript_ID %in% tx_to_keep)

# save FSM and novel matrix to import in seurat
FSM_tx <- DR_data %>%
  dplyr::filter(n_UMI > 0 & structural_category == "full-splice_match") %>%
  pull(transcript_ID)

novel_tx <- DR_data %>%
  dplyr::filter(n_UMI > 0 & structural_category %in% c("novel_in_catalog", "novel_not_in_catalog")) %>%
  pull(transcript_ID)

## count mat with only FSM
FSM_count_mat <-
  umi_filt_count_mat %>%
  dplyr::filter(transcript_ID %in% FSM_tx)
## count mat with FSM, NIC, NNIC
FSM_novel_count_mat <-
  umi_filt_count_mat %>%
  dplyr::filter(transcript_ID %in% c(FSM_tx, novel_tx))
