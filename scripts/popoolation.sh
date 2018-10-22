# need separate pileup for each
cd ~/reciprocal_t/data/aligned

samtools mpileup -Q 20 -B --max-depth 3000 --skip-indels -f /data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.fasta AAAA_F1_REP1.bam > ~/reciprocal_t/analysis/AAAA_F1_REP1.mpileup

samtools mpileup -Q 20 -B --max-depth 3000 --skip-indels -f /data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.fasta AAAA_F1_REP2.bam > ~/reciprocal_t/analysis/AAAA_F1_REP2.mpileup

samtools mpileup -Q 20 -B --max-depth 3000 --skip-indels -f /data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.fasta AAAA_F1_REP3.bam > ~/reciprocal_t/analysis/AAAA_F1_REP3.mpileup

samtools mpileup -Q 20 -B --max-depth 3000 --skip-indels -f /data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.fasta AAAA_F1_REP4.bam > ~/reciprocal_t/analysis/AAAA_F1_REP4.mpileup

samtools mpileup -Q 20 -B --max-depth 3000 --skip-indels -f /data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.fasta HHHH_F1_REP1.bam > ~/reciprocal_t/analysis/HHHH_F1_REP1.mpileup

samtools mpileup -Q 20 -B --max-depth 3000 --skip-indels -f /data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.fasta HHHH_F1_REP2.bam > ~/reciprocal_t/analysis/HHHH_F1_REP2.mpileup

samtools mpileup -Q 20 -B --max-depth 3000 --skip-indels -f /data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.fasta HHHH_F1_REP3.bam > ~/reciprocal_t/analysis/HHHH_F1_REP3.mpileup

samtools mpileup -Q 20 -B --max-depth 3000 --skip-indels -f /data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.fasta HHHH_F1_REP4.bam > ~/reciprocal_t/analysis/HHHH_F1_REP4.mpileup


#### calculate pi

perl ~/bin/popoolation_1.2.2/Variance-at-position.pl --pool-size 20 --min-qual 20 --min-coverage 20 --min-count 2 --max-coverage 500 --min-covered-fraction 0.4 --pileup ~/reciprocal_t/analysis/AAAA_F1_REP1.mpileup --gtf /data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript_NOSPACE.gtf --output AAAA_F1_REP1.genes.pi --measure pi

perl ~/bin/popoolation_1.2.2/Variance-at-position.pl --pool-size 20 --min-qual 20 --min-coverage 20 --min-count 2 --max-coverage 500 --min-covered-fraction 0.4 --pileup ~/reciprocal_t/analysis/AAAA_F1_REP2.mpileup --gtf /data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript_NOSPACE.gtf --output AAAA_F1_REP2.genes.pi --measure pi

perl ~/bin/popoolation_1.2.2/Variance-at-position.pl --pool-size 20 --min-qual 20 --min-coverage 20 --min-count 2 --max-coverage 500 --min-covered-fraction 0.4 --pileup ~/reciprocal_t/analysis/AAAA_F1_REP3.mpileup --gtf /data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript_NOSPACE.gtf --output AAAA_F1_REP3.genes.pi --measure pi

perl ~/bin/popoolation_1.2.2/Variance-at-position.pl --pool-size 20 --min-qual 20 --min-coverage 20 --min-count 2 --max-coverage 500 --min-covered-fraction 0.4 --pileup ~/reciprocal_t/analysis/AAAA_F1_REP4.mpileup --gtf /data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript_NOSPACE.gtf --output AAAA_F1_REP4.genes.pi --measure pi

perl ~/bin/popoolation_1.2.2/Variance-at-position.pl --pool-size 20 --min-qual 20 --min-coverage 20 --min-count 2 --max-coverage 500 --min-covered-fraction 0.4 --pileup ~/reciprocal_t/analysis/HHHH_F1_REP1.mpileup --gtf /data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript_NOSPACE.gtf --output HHHH_F1_REP1.genes.pi --measure pi

perl ~/bin/popoolation_1.2.2/Variance-at-position.pl --pool-size 20 --min-qual 20 --min-coverage 20 --min-count 2 --max-coverage 500 --min-covered-fraction 0.4 --pileup ~/reciprocal_t/analysis/HHHH_F1_REP2.mpileup --gtf /data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript_NOSPACE.gtf --output HHHH_F1_REP2.genes.pi --measure pi

perl ~/bin/popoolation_1.2.2/Variance-at-position.pl --pool-size 20 --min-qual 20 --min-coverage 20 --min-count 2 --max-coverage 500 --min-covered-fraction 0.4 --pileup ~/reciprocal_t/analysis/HHHH_F1_REP3.mpileup --gtf /data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript_NOSPACE.gtf --output HHHH_F1_REP3.genes.pi --measure pi

perl ~/bin/popoolation_1.2.2/Variance-at-position.pl --pool-size 20 --min-qual 20 --min-coverage 20 --min-count 2 --max-coverage 500 --min-covered-fraction 0.4 --pileup ~/reciprocal_t/analysis/HHHH_F1_REP4.mpileup --gtf /data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript_NOSPACE.gtf --output HHHH_F1_REP4.genes.pi --measure pi




#### calculate pi in 1 kb windows

perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --pool-size 40 --min-qual 20 --min-coverage 20 --min-count 2 --max-coverage 500 --min-covered-fraction 0.4 --input ~/reciprocal_t/analysis/AAAA_F1_REP1.mpileup --output AAAA_F1_REP1.1kb.pi --measure pi --window-size 1000 --step-size 500

perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --pool-size 40 --min-qual 20 --min-coverage 20 --min-count 2 --max-coverage 500 --min-covered-fraction 0.4 --input ~/reciprocal_t/analysis/AAAA_F1_REP2.mpileup --output AAAA_F1_REP2.1kb.pi --measure pi --window-size 1000 --step-size 500

perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --pool-size 40 --min-qual 20 --min-coverage 20 --min-count 2 --max-coverage 500 --min-covered-fraction 0.4 --input ~/reciprocal_t/analysis/AAAA_F1_REP3.mpileup --output AAAA_F1_REP3.1kb.pi --measure pi --window-size 1000 --step-size 500

perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --pool-size 40 --min-qual 20 --min-coverage 20 --min-count 2 --max-coverage 500 --min-covered-fraction 0.4 --input ~/reciprocal_t/analysis/AAAA_F1_REP4.mpileup --output AAAA_F1_REP4.1kb.pi --measure pi --window-size 1000 --step-size 500

perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --pool-size 40 --min-qual 20 --min-coverage 20 --min-count 2 --max-coverage 500 --min-covered-fraction 0.4 --input ~/reciprocal_t/analysis/HHHH_F1_REP1.mpileup --output HHHH_F1_REP1.1kb.pi --measure pi --window-size 1000 --step-size 500

perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --pool-size 40 --min-qual 20 --min-coverage 20 --min-count 2 --max-coverage 500 --min-covered-fraction 0.4 --input ~/reciprocal_t/analysis/HHHH_F1_REP2.mpileup --output HHHH_F1_REP2.1kb.pi --measure pi --window-size 1000 --step-size 500

perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --pool-size 40 --min-qual 20 --min-coverage 20 --min-count 2 --max-coverage 500 --min-covered-fraction 0.4 --input ~/reciprocal_t/analysis/HHHH_F1_REP3.mpileup --output HHHH_F1_REP3.1kb.pi --measure pi --window-size 1000 --step-size 500

perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --pool-size 40 --min-qual 20 --min-coverage 20 --min-count 2 --max-coverage 500 --min-covered-fraction 0.4 --input ~/reciprocal_t/analysis/HHHH_F1_REP4.mpileup --output HHHH_F1_REP4.1kb.pi --measure pi --window-size 1000 --step-size 500

