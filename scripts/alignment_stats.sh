cd ~/reciprocal_t/data/aligned

for i in $(ls *.bam | cut -f -1 -d "." | uniq )

do {
ONE=$(samtools view -F 4 ${i}.bam | wc -l) #mapped reads
TWO=$(samtools view -F 4 -q 20 ${i}.bam | wc -l) #mapped reads with quality of 20
THREE=$(samtools view -f 1024 ${i}.bam | wc -l) #pcr duplicates
FOUR=$(samtools view  ${i}.bam | wc -l) #total

# now coverage
samtools depth ${i}.bam | head -n 800000  > ${i}.coverage

shuf ${i}.coverage > ${i}.coverage1
mv ${i}.coverage1 ${i}.coverage

AVG=$(cat ${i}.coverage  | awk '{sum+=$3; sumsq+=$3*$3} END { print sum/NR}')

STD=$(cat ${i}.coverage  | awk '{sum+=$3; sumsq+=$3*$3} END { print sqrt(sumsq/NR - (sum/NR)**2)}')

rm ${i}.coverage

echo ${i},$ONE,$TWO,$THREE,$FOUR,$AVG,$STD

 } >> ~/reciprocal_t/analysis/count.aligned.txt

done
