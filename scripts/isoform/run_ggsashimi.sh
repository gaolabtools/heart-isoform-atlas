#!/bin/bash

# OUTFILE="plots/fig3_TNNI3"
# POS="chr19:55151800-55157800"
# TXID="TNNI3-201,TNNI3-206,TNNI3-208"
# MINR=150
# BAM="normal_bam_lists/normal_template.tsv"

# OUTFILE="plots/fig3_ACTG1"
# POS="chr17:81509900-81512900"
# TXID="ACTG1-205,ACTG1-206,ACTG1-212,ACTG1-204,ACTG1-219,ACTG1-214"
# BAM="normal_bam_lists/normal_template.tsv"
# MINR=10

GTF="GRCh38-2024-A.gtf"

eval "$(conda shell.bash hook)"
conda activate /home/tpg8911/miniconda3/envs/ggsashimi

python /home/tools/ggsashimi/ggsashimi.py \
-b $BAM \
-c $POS \
-M $MINR \
-g $GTF \
--txIDs $TXID \
--ann-height 2.25 \
--aggr mean_j \
--height 1.75 \
--width 10 \
--base-size 24 \
-P palette_onecolor.txt \
-C 3 \
-O 3 \
--alpha 0.2 \
-o $OUTFILE