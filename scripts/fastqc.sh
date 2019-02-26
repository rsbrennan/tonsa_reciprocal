output=~/reciprocal_t/analysis/fastqc/

for i in $(ls /data/copepods/reciprocal_tonsa/raw_data/hwftp.novogene.com/C202SC18061792/raw_data )

do {
    cd /data/copepods/reciprocal_tonsa/raw_data/hwftp.novogene.com/C202SC18061792/raw_data${i}

    ~/bin/FastQC/fastqc *fq.gz -t 4 -f fastq --noextract -o ${output}

    echo ${i} done

  }
done >&2 | tee -a ~/reciprocal_t/log_out/fastqc.error.out

# aggregate fastqc across files:

multiqc ~/reciprocal_t/analysis/fastqc/
