# summary of A. tonsa F21 analysis

These results are from RNAseq data from four replicates of each AAAA and HHHH lines. 

## methods

Variants called with [varscan](http://dkoboldt.github.io/varscan/)
  
  1. All samps have at least 40x coverage
  2. where we call at least 3/8 lines variable.
  3. max coverage per line: Median coverage is 128, remove those > 3*128=385x. But looking across all lines, so across all 8 lines where coverage is > 3078x
  4. This leaves 363,200 variable sites with new median coverage of 115x
  5. 11,026 genes remain (out of 33,313). If we only include genes with at least 10 variants, 7,905 remain. These 7,905 genes are what I'm using from this point forward. 

## heterozygosity

Calculated this a number of ways. All results tell the same story. Only showing one approach for simplicity/clarity.

Looking at heterozygosity across each gene. 

calculated as 	`Hp <- (2 * sum(numMajor) * sum(numMinor)) / ((sum(numMajor) + sum(numMinor))^2)`

where `sum(numMajor)` and `sum(numMinor))` are the replicate specific sums of the count of the reference and alternate allele. These are counted at all snps in the gene. 

Across the genome, there is no reduction in Hp in the HHHH line. 

| Line    |  mean heterozygosity  |
|:------:| :-----:|
| AAAA_F1_REP1 | 0.302 |
| AAAA_F1_REP2 | 0.303 |
| AAAA_F1_REP3 | 0.304 |
| AAAA_F1_REP4 | 0.302 |
| HHHH_F1_REP1 | 0.302 |
| HHHH_F1_REP2 | 0.302 |
| HHHH_F1_REP3 | 0.304 |
| HHHH_F1_REP4 | 0.303 |


<img src="https://raw.githubusercontent.com/rsbrennan/tonsa_f20/master/figures/het_boxplot.png" width="500">


<img src="https://raw.githubusercontent.com/rsbrennan/tonsa_f20/master/figures/heterozygosity.png" width="500">


The low end of quantiles also do not differ:

<img src="https://raw.githubusercontent.com/rsbrennan/tonsa_f20/master/figures/het.01_quantile.png" width="500">

<img src="https://raw.githubusercontent.com/rsbrennan/tonsa_f20/master/figures/het.05_quantile.png" width="500">


Looking at the lower 0.05 quantile of heterozygosity for each line. So lowest 395 genes

The follwoing table is the number of genes that overlap between any 2, 3, or 4 lines in each line

| Shared between N reps | AAAA | HHHH|
|:------:| :-----:|:-----:|
| **Two**   | 415 | 432 |
| **Three** | 253 | 271 |
| **Four**  | 138 | 146 |


87 are consistent across all lines- both AAAA and HHHH.


**There is no systematic reduction in diversity in the HHHH line**








