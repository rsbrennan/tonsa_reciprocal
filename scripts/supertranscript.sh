
cd /data/atonsa/Atonsa_gen_trans_agp_gff/

# collapse transcripts to genes
~/bin/trinityrnaseq-Trinity-v2.8.0/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py \
       --trinity_fasta trinity.trimmomatic.above500.noPhiX.fasta \
       --out_prefix atonsa_super_transcript

~/bin/bwa/bwa index atonsa_super_transcript.fasta
