#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate /home/miniconda3/envs/SQANTI3.env

GTF="OUT.extended_annotation.gtf"
ANNOTATION="ref/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"
GENOME="ref/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
CAGE="tools/SQANTI3-5.2.2/data/ref_TSS_annotation/refTSS_v4.1_human_coordinate.hg38.bed"

python /home/tools/SQANTI3-5.2.2/sqanti3_qc.py $GTF $ANNOTATION $GENOME --CAGE_peak $CAGE -d isoquantv342_refTSSv41 --cpus 8 --report both