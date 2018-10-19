#!/bin/bash -l

my_bwa=~/bin/bwa/bwa
my_samtools=~/bin/samtools-1.6/samtools
bwagenind=/data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.fasta
my_samblstr=~/bin/samblaster/samblaster

# index reference
#$my_bwa index oyster_assembly.fa

cd ~/reciprocal_t/data/trimmed/

for sample in `ls ~/reciprocal_t/data/trimmed | grep '.fq.gz' | cut -f 1-3 -d "_"| uniq | grep -v 'AAAA\|HHHH'
`
do

    echo "starting sample ${sample}"
    #starting with only name of rep. need to pull out files

    lane1_1=$(ls ~/reciprocal_t/data/trimmed | grep ${sample} | grep '_L1_1')
    lane1_2=$(ls ~/reciprocal_t/data/trimmed | grep ${sample} | grep '_L1_2')
    lane2_1=$(ls ~/reciprocal_t/data/trimmed | grep ${sample} | grep '_L2_1')
    lane2_2=$(ls ~/reciprocal_t/data/trimmed | grep ${sample} | grep '_L2_2')

    # align lane 1
    echo "starting lane 1 for ${sample}"

    lib1=$(echo $lane1_1 | cut -f 1-3,6 -d "_")

    rg=$(echo \@RG\\tID:$sample\\tPL:Illumina\\tPU:x\\tLB:$lib1\\tSM:$sample)

    $my_bwa mem -t 4 -R $rg $bwagenind $lane1_1 $lane1_2 | \
    $my_samblstr |\
    $my_samtools view -h -u - | \
    $my_samtools sort - -O bam -o ~/reciprocal_t/data/aligned/${lib1}.bam

    # align lane 2
    echo "starting lane 2 for ${sample}"

    lib2=$(echo $lane2_1 | cut -f 1-3,6 -d "_")

    rg=$(echo \@RG\\tID:$sample\\tPL:Illumina\\tPU:x\\tLB:$lib2\\tSM:$sample)

    $my_bwa mem -t 4 -R $rg $bwagenind $lane2_1 $lane2_2 | \
    $my_samblstr |\
    $my_samtools view -h -u - | \
    $my_samtools sort - -O bam -o ~/reciprocal_t/data/aligned/${lib2}.bam

    # merge the two bam files

    ~/bin/bamtools/bamtools merge -in ~/reciprocal_t/data/aligned/${lib1}.bam -in ~/reciprocal_t/data/aligned/${lib2}.bam | $my_samtools sort - -O bam -o ~/reciprocal_t/data/aligned/${sample}.bam

done

