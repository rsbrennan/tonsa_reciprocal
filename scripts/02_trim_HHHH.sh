cd ~/reciprocal_t/data/trimmed

#cp /data/programs/Trimmomatic-0.36/adapters/TruSeq3-PE.fa .

for i in $(ls -d /data/reciprocal_t/raw_data/hwftp.novogene.com/C202SC18061792/raw_data/HHHH* | cut -f 8 -d "/" )

do {

ls /data/reciprocal_t/raw_data/hwftp.novogene.com/C202SC18061792/raw_data/${i}/*.fq.gz | grep '_L1_1.fq.gz'


    read1=$(ls /data/reciprocal_t/raw_data/hwftp.novogene.com/C202SC18061792/raw_data/${i}/*.fq.gz | grep '_L1_1.fq.gz')
    read2=$(ls /data/reciprocal_t/raw_data/hwftp.novogene.com/C202SC18061792/raw_data/${i}/*.fq.gz | grep '_L1_2.fq.gz')
    base1=$(ls /data/reciprocal_t/raw_data/hwftp.novogene.com/C202SC18061792/raw_data/${i}/*.fq.gz | grep '_L1_1.fq.gz' | cut -f 9- -d "/" | cut -f 1 -d ".")
    base2=$(ls /data/reciprocal_t/raw_data/hwftp.novogene.com/C202SC18061792/raw_data/${i}/*.fq.gz | grep '_L1_2.fq.gz' | cut -f 9- -d "/" | cut -f 1 -d ".")


  java -jar /data/programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
   -threads 2 \
    ${read1} ${read2} \
    ${base1}.qc.fq.gz s1_sed \
    ${base2}.qc.fq.gz s2_sed \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:TRUE \
    LEADING:2 TRAILING:2 \
    SLIDINGWINDOW:4:2 \
    MINLEN:31

    echo $i lane 1 done

    ### run lane 2

    ls /data/reciprocal_t/raw_data/hwftp.novogene.com/C202SC18061792/raw_data/${i}/*.fq.gz | grep '_L2_1.fq.gz'


    read1=$(ls /data/reciprocal_t/raw_data/hwftp.novogene.com/C202SC18061792/raw_data/${i}/*.fq.gz | grep '_L2_1.fq.gz')
    read2=$(ls /data/reciprocal_t/raw_data/hwftp.novogene.com/C202SC18061792/raw_data/${i}/*.fq.gz | grep '_L2_2.fq.gz')
    base1=$(ls /data/reciprocal_t/raw_data/hwftp.novogene.com/C202SC18061792/raw_data/${i}/*.fq.gz | grep '_L2_1.fq.gz' | cut -f 9- -d "/" | cut -f 1 -d ".")
    base2=$(ls /data/reciprocal_t/raw_data/hwftp.novogene.com/C202SC18061792/raw_data/${i}/*.fq.gz | grep '_L2_2.fq.gz' | cut -f 9- -d "/" | cut -f 1 -d ".")


  java -jar /data/programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
   -threads 2 \
    ${read1} ${read2} \
    ${base1}.qc.fq.gz s1_sed \
    ${base2}.qc.fq.gz s2_sed \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:TRUE \
    LEADING:2 TRAILING:2 \
    SLIDINGWINDOW:4:2 \
    MINLEN:31

    echo $i lane 2 done

  }

done
