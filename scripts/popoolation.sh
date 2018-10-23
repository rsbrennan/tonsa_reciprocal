# need separate pileup for each

for i in $(ls ~/reciprocal_t/data/aligned | cut -f 1 -d '.'| grep 'AAHH'); do

    echo "Status: starting $i"

    samtools mpileup -Q 20 -B --max-depth 2000 --skip-indels \
    -f /data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.fasta \
    ~/reciprocal_t/data/aligned/${i}.bam \
    > ~/reciprocal_t/analysis/popoolation/${i}.mpileup

    echo "Status: $i pileup done; starting popoolation"

    perl ~/bin/popoolation_1.2.2/Variance-at-position.pl --pool-size 20 --min-qual 20 \
    --min-coverage 20 --min-count 2 --max-coverage 500 --min-covered-fraction 0.4 \
    --pileup ~/reciprocal_t/analysis/popoolation/${i}.mpileup \
    --gtf /data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript_NOSPACE.gtf \
    --output ~/reciprocal_t/analysis/popoolation/${i}.genes.pi --measure pi

    rm ~/reciprocal_t/analysis/popoolation/${i}.mpileup

    echo "Status: $i popoolation done, pileup removed"

done
