
# there are quotes in these files, remove them. Could do it earlier, but.. I haven't.

for i in snp_F1 snp_F2 snp_F3 dge_F1 dge_F2 dge_F3; do

  echo ${i}_GOterms.out

  cat ~/reciprocal_t/analysis/GO_enrich/${i}_GOterms.out | sed 's/"//g' > ~/reciprocal_t/analysis/GO_enrich/${i}_GOterms.corrected.out

done

# make dge files comme delim. again, should have done this when generating the files... I didn't.

for i in F1 F2 F3; do

  cat ~/reciprocal_t/analysis/GO_enrich/dge_contrib_${i}.txt | cut -f 1 > ~/reciprocal_t/analysis/GO_enrich/dge_contrib_${i}.genes
  cat ~/reciprocal_t/analysis/GO_enrich/dge_contrib_${i}.txt  | cut -f 2 | paste ~/reciprocal_t/analysis/GO_enrich/dge_contrib_${i}.genes - -d ","  > ~/reciprocal_t/analysis/GO_enrich/dge_contrib_${i}.go.in

done

