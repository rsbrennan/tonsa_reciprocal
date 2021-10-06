#!/bin/bash -l

cd ~/reciprocal_t/analysis/SNPvalidation/merged/

samtools mpileup -Q 20 -B --max-depth 3000 --skip-indels -f /data/copepods/tonsa_transcriptome/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.fasta AA_F25_Rep1.bam  AA_F25_Rep3.bam  HH_F25_Rep1.bam  HH_F25_Rep3.bam AA_F25_Rep2.bam  AA_F25_Rep4.bam  HH_F25_Rep2.bam  HH_F25_Rep4.bam | java -jar ~/bin/VarScan.v2.3.9.jar mpileup2snp --mpileup 1 --min-coverage 30 --min-reads 10 --min-avg-qual 20 --min-var-freq 0.01 --variants --p-value 0.1 > ~/reciprocal_t/analysis/SNPvalidation/snp_validation_out

