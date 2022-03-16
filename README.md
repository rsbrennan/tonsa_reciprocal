#Acartia tonsa reciprocal transplant

This repository hold scripts for the analysis of a reciprocal transplant of Acartia tonsa following a ~20 generation selection experiment.  

These RNAseq data were sequenced on 2 lanes of a Novaseq by Novogene. Data were received by us on 2018-07-24. 

## Data availability

Raw sequence data is available at NCBI BioProject PRJNA555881

Gene expression data are included in this repository in `DGE_data/`, both normalized (by librarysize) and not normalized.

All other summaries etc. are included as supplemental files with the manuscript.

Please get in touch about any data issues.

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
- Generate the correlation scatter plot for DGE and AF for the supplemental: `13_scatter_plot.R`
- GO analysis
  - formatted output for DGE is from `13_scatter_plot.R`: "~/reciprocal_t/analysis/GO_enrich/dge_f1.txt"
  - `14b_go_assign_snps.py`, `14b_go_assign_dge.py` produces `dge_F1_GOterms.out` `snp_GOterms.out`, for each generation.
  - `14c_GO_format.sh`there are some weird formatting issues (some quotes?) that were just easier to fix with bash
  - `14d_GO_Enrich.R` and `14d_GO_Enrich_SNP.R` run the actual GO enrichment for expression and snps, respectively
  - `14e_go_enrich_deltapi.md` run go enrichment for the change in pi

## Figures:

- Fig. 2: combination of `Fig_pca.R` and `14d_GO_Enrich.R`  
- Fig. 3: `plasticity_analysis.R`   
- Fig. 4: `11_dapc.R`    
- Fig. 5: `10b_popoolation_pi.R`
- Fig. 6: `15_Phenotype_data_analysis.R`

Supplemental:
- Fig. S1: `Fig_S1.R`  
- Fig. S2: `indiv_gene_fig.md`  
- Fig. S3: `indiv_gene_fig.md`  
- Fig. S4: `13_scatter_plot.R`  
- Fig. S5: `snp_validation.md`  
- Fig. S6: `13_scatter_plot.R`  
- Fig. S7: `09_filter_variants.R`  
- Fig. S8: `venn.R`  
- Fig. S9: `14d_GO_Enrich.R`  
- Fig. S10: `14d_GO_Enrich.R`  
