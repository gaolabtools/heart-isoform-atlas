#! /bin/bash

P_DIR="/home/tools/scNanoGPS_v0.15/"
FASTQ="/data/nanopore/fastq_pass/"
ISOQUANT="/home/miniconda3/envs/isoquant_v3.4.2/bin/isoquant.py"
REF_GENOME="/home/tools/genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
IND_GENOME="/home/tools/genome/refdata-gex-GRCh38-2024-A/fasta/refdata-gex-GRCh38-2024-A.mmi"
GENOME_ANNOTATION="/home/tools/genome/refdata-gex-GRCh38-2024-A/genes/genes.gtf"
rRNA_HB_BED="/home/tools/gencode/rRNA_HB.gencode_v44.gene.bed"
ncores=20
ANNOVAR="/home/tools/annovar"
ANNOVAR_DB="/home/tools/annovar/hg38db/"
ANNOVAR_GV="hg38"
ANNOVAR_PROTOCOL="refGene,cytoBand,gnomad30_genome,avsnp150,dbnsfp42c,cosmic96_coding,cosmic96_noncoding"
ANNOVAR_OP="gx,r,f,f,f,f,f"
ANNOVAR_XREF="/home/tools/annovar/hg38db/omim/gene_xref.txt"

python3 $P_DIR/other_utils/read_length_profiler.py -i $FASTQ &> run_read_length_profiler.log.txt &
python3 $P_DIR/scanner.py -i $FASTQ -t $ncores &> run_scanner.log.txt
python3 $P_DIR/assigner.py -t $ncores &> run_assigner.log.txt
python3 $P_DIR/curator.py -t $ncores --ref_genome $REF_GENOME --idx_genome $IND_GENOME --exc_bed $rRNA_HB_BED &> run_curator.log.txt
python3 $P_DIR/reporter_expression.py --gtf $GENOME_ANNOTATION -t $ncores &> run_reporter_expression.log.txt
python3 $P_DIR/reporter_isoform.isoquant.py --ref_genome $REF_GENOME --gtf $GENOME_ANNOTATION -t $ncores --isoquant $ISOQUANT &> run_reporter_isoform.log.txt
python3 $P_DIR/reporter_SNV.py --ref_genome $REF_GENOME -t $ncores --annovar $ANNOVAR --annovar_db $ANNOVAR_DB --annovar_gv $ANNOVAR_GV --annovar_protocol $ANNOVAR_PROTOCOL --annovar_operation $ANNOVAR_OP --annovar_xref $ANNOVAR_XREF &> run_reporter_SNV.log.txt
python3 $P_DIR/other_utils/parse_annovar_column.py -i scNanoGPS_res/annovar.hg38_multianno.vcf > scNanoGPS_res/annovar.hg38_multianno.tsv
python3 $P_DIR/reporter_summary.py --ref_genome $REF_GENOME --gtf $GENOME_ANNOTATION --qualimap_param "--java-mem-size=300G" &> run_reporter_summary.log.txt

