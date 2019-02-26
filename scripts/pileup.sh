
# pipe it
#samtools mpileup -Q 20 -B --max-depth 3000 --skip-indels -f /data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.fasta AAAA_F1_REP1.bam AAAA_F1_REP2.bam AAAA_F1_REP3.bam AAAA_F1_REP4.bam HHHH_F1_REP1.bam HHHH_F1_REP2.bam HHHH_F1_REP3.bam HHHH_F1_REP4.bam | java -jar ~/bin/VarScan.v2.3.9.jar mpileup2snp --mpileup 1 --min-coverage 40 --min-reads 2 --min-avg-qual 20 --min-var-freq 0.025 --variants --p-value 0.1 > ~/reciprocal_t/analysis/snp_out

# all reps

cd ~/reciprocal_t/data/aligned

samtools mpileup -Q 20 -B --max-depth 3000 --skip-indels -f /data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.fasta AAAA_F1_REP1.bam AAAA_F1_REP2.bam AAAA_F1_REP3.bam AAAA_F1_REP4.bam AAAA_F2_REP1.bam AAAA_F2_REP2.bam AAAA_F2_REP3.bam AAAA_F2_REP4.bam AAAA_F3_REP1.bam AAAA_F3_REP2.bam AAAA_F3_REP3.bam AAAA_F3_REP4.bam AAHH_F1_REP1.bam AAHH_F1_REP2.bam AAHH_F1_REP3.bam AAHH_F1_REP4.bam AAHH_F2_REP1.bam AAHH_F2_REP2.bam AAHH_F2_REP3.bam AAHH_F2_REP4.bam AAHH_F3_REP1.bam AAHH_F3_REP2.bam AAHH_F3_REP3.bam AAHH_F3_REP4.bam HHAA_F1_REP1.bam HHAA_F1_REP2.bam HHAA_F1_REP3.bam HHAA_F1_REP4.bam HHAA_F2_REP1.bam HHAA_F2_REP2.bam HHAA_F2_REP3.bam HHAA_F2_REP4.bam HHAA_F3_REP1.bam HHAA_F3_REP2.bam HHAA_F3_REP3.bam HHAA_F3_REP4.bam HHHH_F1_REP1.bam HHHH_F1_REP2.bam HHHH_F1_REP3.bam HHHH_F1_REP4.bam HHHH_F2_REP1.bam HHHH_F2_REP2.bam HHHH_F2_REP3.bam HHHH_F2_REP4.bam HHHH_F3_REP1.bam HHHH_F3_REP2.bam HHHH_F3_REP3.bam HHHH_F3_REP4.bam | java -jar ~/bin/VarScan.v2.3.9.jar mpileup2snp --mpileup 1 --min-coverage 30 --min-reads 10 --min-avg-qual 20 --min-var-freq 0.01 --variants --p-value 0.1 > ~/reciprocal_t/analysis/snp_all_out

