#!/bin/bash -l

java -jar ~/bin/VarScan.v2.3.9.jar mpileup2snp ~/reciprocal_t/analysis/merged.pileup --min-coverage 30 --min-reads 2 --min-avg-qual 20 --min-var-freq 0.01 --p-value 0.1 > ~/reciprocal_t/analysis/snp_out

