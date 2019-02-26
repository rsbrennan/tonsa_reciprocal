# will use the quasi mapping approach- it is faster than the original. It is actually unnecessary to specify as it is the default
# here, k is the minimum acceptable length for a valid match. This works well with 75 bp reads or longer

salmon index -t /data/atonsa/Atonsa_gen_trans_agp_gff/trinity.trimmomatic.above500.noPhiX.fasta.gz \
    -i /data/atonsa/Atonsa_gen_trans_agp_gff/tonsa_index --type quasi -k 31
