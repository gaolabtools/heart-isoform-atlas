#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate /home/miniconda3/envs/isoquant_v3.4.2

# Run isoquant pipeline
isoquant.py --reference ref/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
--genedb ref/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz \
--bam_list bam_list.txt \
--read_group tag:RG \
--data_type nanopore \
--complete_genedb \
--sqanti_output \
--threads 12 \
-o results/v3.4.2