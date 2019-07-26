# Acartia tonsa reciprocal transplant

This repository hold scripts for the analysis of a reciprocal transplant of Acartia tonsa following a ~20 generation selection experiment.  

These RNAseq data were sequenced on 2 lanes of a Novaseq by Novogene. Data were received by us on 2018-07-24. 

## Data availability

Raw sequence data is available at NCBI BioProject PRJNA555881

All other summaries etc. are included as supplemental files with the manuscript.

## Scripts

Below are scripts to run the full analysis for the manuscript. A short description accompanies each. 

### Data processing

- Check data quality: `01_fastqc.sh`  
- Trim fastq files: `02_trim_AAAA.sh`; `02_trim_HHHH.sh`;`02_trim_HHAA.sh`;`02_trim_AAHH.sh`
- Re-run fastqc: `03_fastqc_posttrim.sh`

### Aligning, read counts, variant calling

- create index for salmon: `04_salmon_index.sh`
- Quantify each sample: `05_salmon.sh` 
- Generate supertranscript reference to align: `06_supertranscript.sh`
- Align to the supertranscript: `07_align.sh`
- call variants using varscan: `08_varscan_all.sh`
  - note that this is very liberal in the calls. Need to filter.
- Filter variants from varscan: `09_filter_variants.R`

### Analysis

## Figures:

Fig. 1:  
Fig. 2:  
Fig. 3:
Fig. 4:
Fig. S1: 


