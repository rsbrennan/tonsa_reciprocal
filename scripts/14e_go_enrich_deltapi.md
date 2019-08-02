
```python
import pandas as pd
import numpy as np
import gzip
import csv
import re
from collections import OrderedDict


print("Starting SNP go assignments")

# assign GO terms
filepath = '/data/atonsa/Atonsa_gen_trans_agp_gff/blast2go_gff_export_20170613_0815.gff.gz'

####
## SNPs pi_aahh.
####

#make empty array
num_lines = sum(1 for line in open('/users/r/b/rbrennan/reciprocal_t/analysis/GO_enrich/pi_aahh.txt'))-1

go_df = np.empty(shape=(num_lines,2), dtype = object)

ict=0 # start counter

with open('/users/r/b/rbrennan/reciprocal_t/analysis/GO_enrich/pi_aahh.txt') as master_file:
        head1 = next(master_file)
        for idx, line in enumerate(master_file):
            tmp_gene = line.split(",")[0]
            out_go = []
            # match GO terms
            match = []
            with gzip.open(filepath, mode="rt") as file:
                for gene_line in file:
                    if len(re.findall(tmp_gene, gene_line)) > 0:
                        gene_spl = gene_line.split("\n")[0].split("\t")
                        go_spl = [i for i in gene_spl[8].split(";") if i.startswith('Ontology_id')]
                        if len(go_spl) > 0:
                            if len(go_spl[0].split()) > 1:
                                go_tmp = go_spl[0].split()[1]
                            if len(out_go) == 0:
                                out_go = go_tmp
                            if len(out_go) > 0:
                                out_go = out_go + "," + go_tmp
            go_df[idx,0] = tmp_gene
            if len(out_go) > 0:
                out_go = ";".join(set(out_go.split(","))) # remove duplicates, join together with string.
                go_df[idx,1] = out_go
                #go_df[idx,1] = out_go.replace('"',"").replace(',',";")
            if len(out_go) == 0:
                go_df[idx,1] = "unknown"
            ict=ict+1
            if ict % 100 == 0: print(ict)

print("F1 done")

np.savetxt('/users/r/b/rbrennan/reciprocal_t/analysis/GO_enrich/pi_aahh_GOterms.out', go_df,fmt='%s', delimiter='\t')

print("F1 saved")




####
## SNPs pi_hhaa.
####

#make empty array
num_lines = sum(1 for line in open('/users/r/b/rbrennan/reciprocal_t/analysis/GO_enrich/pi_hhaa.txt'))-1

go_df = np.empty(shape=(num_lines,2), dtype = object)

ict=0 # start counter

with open('/users/r/b/rbrennan/reciprocal_t/analysis/GO_enrich/pi_hhaa.txt') as master_file:
        head1 = next(master_file)
        for idx, line in enumerate(master_file):
            tmp_gene = line.split(",")[0]
            out_go = []
            # match GO terms
            match = []
            with gzip.open(filepath, mode="rt") as file:
                for gene_line in file:
                    if len(re.findall(tmp_gene, gene_line)) > 0:
                        gene_spl = gene_line.split("\n")[0].split("\t")
                        go_spl = [i for i in gene_spl[8].split(";") if i.startswith('Ontology_id')]
                        if len(go_spl) > 0:
                            if len(go_spl[0].split()) > 1:
                                go_tmp = go_spl[0].split()[1]
                            if len(out_go) == 0:
                                out_go = go_tmp
                            if len(out_go) > 0:
                                out_go = out_go + "," + go_tmp
            go_df[idx,0] = tmp_gene
            if len(out_go) > 0:
                out_go = ";".join(set(out_go.split(","))) # remove duplicates, join together with string.
                go_df[idx,1] = out_go
                #go_df[idx,1] = out_go.replace('"',"").replace(',',";")
            if len(out_go) == 0:
                go_df[idx,1] = "unknown"
            ict=ict+1
            if ict % 100 == 0: print(ict)

print("F1 done")

np.savetxt('/users/r/b/rbrennan/reciprocal_t/analysis/GO_enrich/pi_hhaa_GOterms.out', go_df,fmt='%s', delimiter='\t')

print("F1 saved")

```



# run go enrich

```r

setwd("~/reciprocal_t/analysis/GO_enrich/")
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")

hhaa <- read.csv("pi_hhaa.txt", header=TRUE)

cutoff_hhaa <- quantile(hhaa$mean_pi, 0.1)


input="pi_hhaa.txt" # loadings scores from dapc
goAnnotations="pi_hhaa_GOterms.out" #snps

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
    perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
    largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
    smallest=5,   # a GO category should contain at least this many genes to be considered
    clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
   Alternative="l" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
)


#png("~/reciprocal_t/figures/pi_hhaa_GO.png",width = 7, height = 4, res=300, units="in")

results=gomwuPlot(input,goAnnotations,goDivision,
    absValue=cutoff_hhaa,
    level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories 
                # containing genes exceeding the absValue.
    level2=0.05, # FDR cutoff to print in regular (not italic) font.
    level3=0.01,
    txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
    treeHeight=0.5, # height of the hierarchical clustering tree
#   colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
)

#dev.off()

### hh in aa

aahh <- read.csv("pi_aahh.txt", header=TRUE)

cutoff_aahh <- quantile(aahh$mean_pi, 0.05)

input="pi_aahh.txt" # loadings scores from dapc
goAnnotations="pi_aahh_GOterms.out" #snps

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
    perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
    largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
    smallest=5,   # a GO category should contain at least this many genes to be considered
    clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
   Alternative="l" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
)


pdf("~/reciprocal_t/figures/pi_aahh_GO.pdf",width = 8, height = 5)

results=gomwuPlot(input,goAnnotations,goDivision,
    #absValue=cutoff_aahh,
    level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories 
                # containing genes exceeding the absValue.
    level2=0.05, # FDR cutoff to print in regular (not italic) font.
    level3=0.01,
    txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
    treeHeight=0.5, # height of the hierarchical clustering tree
#   colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
)

dev.off()


```

