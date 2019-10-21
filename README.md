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

- calulate pi and look at the loss: 
  - `10a_popoolation.sh`, `10b_popoolation_pi.R`
- run the DAPC: `11_dapc.R`
- Run the CMH, output allele freq changes: `12_snp_analysis.R`
- Generate the scatter plot for Figure 3: `13_scatter_plot.R`
- GO analysis
  - `14a_loadings_for_GO.R` produces formatted snp output, `from DAPC.R`
  - `14b_go_assign_snps.py`, `14b_go_assign_dge.py` produces `dge_F1_GOterms.out` `snp_F1_GOterms.out`, for each generation.
  - `14c_GO_format.sh`there are some weird formatting issues (some quotes?) that were just easier to fix with bash
  - `14d_GO_MWU.R` run the actual GO enrichment
  - `14e_go_enrich_deltapi.md` run go enrichment for the change in pi
- log-rank survivorship, egg production and fecundity calculations, malthusian parameter calculations: `15_Phenotype_data_analysis.R`

## Figures:

Fig. 1: `Fig1_pca.R`  
Fig. 2: `11_dapc.R`  
Fig. 3: `13_scatter_plot.R`  
Fig. 4: `10b_popoolation_pi.R`  
Fig. S1: `Fig_S1.R`  


