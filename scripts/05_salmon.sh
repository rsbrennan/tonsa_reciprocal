
######
#
# quantify each sample with salmon
#
#######

# I have no idea how long this will take.

# -l A tells salmon that it should automatically determine the library type of the sequencing reads (e.g. stranded vs. unstranded etc.)
# -p 8 says uses 8 threads
# -o indicates the directory and name of output
# -l A means automatically infer library type
# seqbias corrects for random hexamer priming
# gcbias corrects for gcbias, but only when present.


for i in $(ls ~/reciprocal_t/data/trimmed | grep '.fq.gz' | cut -f 1-3 -d "_"| uniq);
do

    echo "starting sample ${i}"
    #starting with only name of rep. need to pull out files

    lane1_1=$(ls ~/reciprocal_t/data/trimmed | grep ${i} | grep '_L1_1')
    lane1_2=$(ls ~/reciprocal_t/data/trimmed | grep ${i} | grep '_L1_2')
    lane2_1=$(ls ~/reciprocal_t/data/trimmed | grep ${i} | grep '_L2_1')
    lane2_2=$(ls ~/reciprocal_t/data/trimmed | grep ${i} | grep '_L2_2')

    salmon quant -i /data/atonsa/Atonsa_gen_trans_agp_gff/tonsa_index \
         -l A \
         -1 ~/reciprocal_t/data/trimmed/${lane1_1} ~/reciprocal_t/data/trimmed/${lane2_1} \
         -2 ~/reciprocal_t/data/trimmed/${lane1_2} ~/reciprocal_t/data/trimmed/${lane2_2} \
         -p 8  \
         --seqBias \
         --gcBias \
         -o ~/reciprocal_t/analysis/salmon/${i}

    echo "sample ${i} done"

done
