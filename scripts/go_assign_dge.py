import pandas as pd
import numpy as np
import gzip
import csv
import re
from collections import OrderedDict

# assign GO terms
filepath = '/data/atonsa/Atonsa_gen_trans_agp_gff/blast2go_gff_export_20170613_0815.gff.gz'


print("Starting DGE go assignments")

####
## dge F1
####

#make empty array
num_lines = sum(1 for line in open('/users/r/b/rbrennan/reciprocal_t/analysis/GO_enrich/dge_contrib_F1.txt'))-1

go_df = np.empty(shape=(num_lines,2), dtype = object)

ict=0 # start counter

with open('/users/r/b/rbrennan/reciprocal_t/analysis/GO_enrich/dge_contrib_F1.txt') as master_file:
        head1 = next(master_file)
        for idx, line in enumerate(master_file):
            tmp_gene = line.split("\t")[0]
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
            if ict % 500 == 0: print(ict)

print("F1 done")

np.savetxt('/users/r/b/rbrennan/reciprocal_t/analysis/GO_enrich/dge_F1_GOterms.out', go_df,fmt='%s', delimiter='\t')

print("F1 saved")

####
## dge F2
####

#make empty array
num_lines = sum(1 for line in open('/users/r/b/rbrennan/reciprocal_t/analysis/GO_enrich/dge_contrib_F2.txt'))-1

go_df = np.empty(shape=(num_lines,2), dtype = object)

ict=0 # start counter

with open('/users/r/b/rbrennan/reciprocal_t/analysis/GO_enrich/dge_contrib_F2.txt') as master_file:
        head1 = next(master_file)
        for idx, line in enumerate(master_file):
            tmp_gene = line.split("\t")[0]
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
            if ict % 500 == 0: print(ict)

print("F2 done")

np.savetxt('/users/r/b/rbrennan/reciprocal_t/analysis/GO_enrich/dge_F2_GOterms.out', go_df,fmt='%s', delimiter='\t')

print("F2 saved")

####
## dge F3
####

#make empty array
num_lines = sum(1 for line in open('/users/r/b/rbrennan/reciprocal_t/analysis/GO_enrich/dge_contrib_F3.txt'))-1

go_df = np.empty(shape=(num_lines,2), dtype = object)

ict=0 # start counter

with open('/users/r/b/rbrennan/reciprocal_t/analysis/GO_enrich/dge_contrib_F3.txt') as master_file:
        head1 = next(master_file)
        for idx, line in enumerate(master_file):
            tmp_gene = line.split("\t")[0]
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
            if ict % 500 == 0: print(ict)

print("F3 done")

np.savetxt('/users/r/b/rbrennan/reciprocal_t/analysis/GO_enrich/dge_F3_GOterms.out', go_df,fmt='%s', delimiter='\t')

print("F3 saved")
