my_bwa=~/bin/bwa/bwa
my_samtools=~/bin/samtools-1.6/samtools
bwagenind=/data/copepods/tonsa_transcriptome/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.fasta
my_samblstr=~/bin/samblaster/samblaster

cd ~/tonsa_genomics/data/trimmed/lane2/

for sample in `ls ~/tonsa_genomics/data/trimmed/lane2 | grep '.fq.gz' | cut -f 1 -d "."| uniq | grep -E 'AA_F25|HH_F25'`

do

    echo "starting sample ${sample}"
    #starting with only name of rep. need to pull out files

    rep_1=$(ls ~/tonsa_genomics/data/trimmed/lane2 | grep ${sample} | grep 'R1')
    rep_2=$(ls ~/tonsa_genomics/data/trimmed/lane2 | grep ${sample} | grep 'R2')

    echo $rep_1
    echo $rep_2

    rg=$(echo \@RG\\tID:$sample\\tPL:Illumina\\tPU:x\\tLB:$sample\\tSM:$sample)

    $my_bwa mem -t 4 -R $rg $bwagenind $rep_1 $rep_2 | \
    $my_samtools view -h -u - | \
    $my_samtools sort - -O bam -o ~/reciprocal_t/analysis/SNPvalidation/${sample}.lane2.bam

done
