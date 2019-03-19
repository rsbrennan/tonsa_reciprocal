# Acartia tonsa reciprocal transplant

This repository hold scripts for the analysis of a reciprocal transplant of Acartia tonsa following a ~20 generation selection experiment. Experiment performed by the Dam lab at UConn. 

These RNAseq data were sequenced on 2 lanes of a Novaseq by Novogene. Data were received by us on 2018-07-24. 

## Scripts

### Salmon

- create index for salmon: `04_salmon_index.sh`
- Quantify each sample: `05_salmon.sh` 

### Call SNPs

- Generate supertranscript reference to align: `06_supertranscript.sh`
- Align to the supertranscript: `07_align.sh`
- call variants using varscan: `08_varscan_all.sh`
  - note that this is very liberal in the calls. Need to filter.
- Filter variants from varscan: `09_filter_variants.R`

### DAPC


### DGE


### SNP CMH


## Figures:

Fig1:  
Fig2:  
Fig3:  
