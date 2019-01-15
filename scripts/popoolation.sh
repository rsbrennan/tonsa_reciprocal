# need separate pileup for each

for i in $(ls ~/reciprocal_t/data/aligned | cut -f 1 -d '.'| grep 'HHAA'); do

    echo "Status: starting $i"

    samtools mpileup -Q 20 -B --max-depth 2000 --skip-indels \
    -f /data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.fasta \
    ~/reciprocal_t/data/aligned/${i}.bam \
    > ~/reciprocal_t/analysis/popoolation/${i}.mpileup

    echo "Status: $i pileup done; starting popoolation"

    perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --pool-size 40 --min-qual 20 \
    --min-coverage 40 --min-count 2 --max-coverage 396 --min-covered-fraction 0.6 \
    --input ~/reciprocal_t/analysis/popoolation/${i}.mpileup \
    --window-size 100 --step-size 100 \
    --output ~/reciprocal_t/analysis/popoolation/${i}.100bp.pi --measure pi

    rm ~/reciprocal_t/analysis/popoolation/${i}.mpileup

    echo "Status: $i popoolation done, pileup removed"

done
