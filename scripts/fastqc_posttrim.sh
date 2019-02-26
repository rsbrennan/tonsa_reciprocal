cd ~/reciprocal_t/data/trimmed/

output=~/reciprocal_t/analysis/fastqc/trim

~/bin/FastQC/fastqc *fq.gz -t 10 -f fastq --noextract -o ${output}

# aggregate fastqc across files:

multiqc ~/reciprocal_t/analysis/fastqc/trim
mv multiqc_* ~/reciprocal_t/analysis/fastqc/trim
