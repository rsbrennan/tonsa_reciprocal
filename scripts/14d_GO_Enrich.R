#### determine "good genes"

########################################
########
######## DGE enrichment
########
########################################

library(dplyr)
library(stringr)

setwd("~/reciprocal_t/analysis/GO_enrich")

df <- read.table("~/reciprocal_t/analysis/AF_change.txt", header=TRUE)
df$gene <- str_split_fixed(df$SNP, ":", 2)[,1]


df_mean <- df %>% rowwise() %>%
  mutate(mean_gen = mean(c(pval_f1, pval_f2, pval_f3), na.rm=TRUE))

df_all <-
  df_mean %>%
  group_by(gene) %>%
  dplyr::summarise(mean_all = min(mean_gen, na.rm=TRUE))

d_f1 <-
  df %>%
  group_by(gene) %>%
  dplyr::summarise(f1_pval = min(pval_f1, na.rm=TRUE))

d_f2 <-
  df %>%
  group_by(gene) %>%
  dplyr::summarise(f2_pval = min(pval_f2, na.rm=TRUE))

d_f3 <-
  df %>%
  group_by(gene) %>%
  dplyr::summarise(f3_pval = min(pval_f3, na.rm=TRUE))


genes <- unique(df$gene)


write.table(data.frame(
    x=df_all$gene, 
    y=-log10(df_all$mean_all)), 
    file="~/reciprocal_t/analysis/GO_enrich/snp_all.txt",
    sep=",", quote=FALSE, col.names=F, row.names=F)

write.table(data.frame(
    x=d_f1$gene, 
    y=-log10(d_f1$f1_pval)), 
    file="~/reciprocal_t/analysis/GO_enrich/snp_f1.txt",
    sep=",", quote=FALSE, col.names=F, row.names=F)

write.table(data.frame(
    x=d_f2$gene, 
    y=-log10(d_f2$f2_pval)), 
    file="~/reciprocal_t/analysis/GO_enrich/snp_f2.txt",
    sep=",", quote=FALSE, col.names=F, row.names=F)

write.table(data.frame(
    x=d_f3$gene, 
    y=-log10(d_f3$f3_pval)), 
    file="~/reciprocal_t/analysis/GO_enrich/snp_f3.txt",
    sep=",", quote=FALSE, col.names=F, row.names=F)

snpF1 <- read.csv("~/reciprocal_t/analysis/GO_enrich/snp_f1.txt", header=FALSE, stringsAsFactors=FALSE)
snpF2 <- read.csv("~/reciprocal_t/analysis/GO_enrich/snp_f2.txt", header=FALSE, stringsAsFactors=FALSE)
snpF3 <- read.csv("~/reciprocal_t/analysis/GO_enrich/snp_f3.txt", header=FALSE, stringsAsFactors=FALSE)
cutoff_snpAll <- quantile(-log10(df_all$mean_all), 0.9, na.rm=T)
cutoff_snpF1 <- quantile(snpF1$V2, 0.9, na.rm=T)
cutoff_snpF2 <- quantile(snpF2$V2, 0.9, na.rm=T)
cutoff_snpF3 <- quantile(snpF3$V2, 0.9, na.rm=T)



# F1

goDatabase="go.obo"
goDivision="BP"
source("gomwu.functions.R")

input="snp_all.txt"
goAnnotations="snp_F1_GOterms.corrected.out"

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
    perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
    largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
    smallest=5,   # a GO category should contain at least this many genes to be considered
    clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
   Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
)

png("~/reciprocal_t/figures/snp_all_GO.png",width = 7, height = 8, res=300, units="in")

results=gomwuPlot(input,goAnnotations,goDivision,
    absValue=cutoff_snpAll,
    level1=0.05, # FDR threshold for plotting. Specify level1=1 to plot all GO categories 
    level2=0.05, # FDR cutoff to print in regular (not italic) font.
    level3=0.01,
    txtsize=1.2,    
    treeHeight=0.5, # height of the hierarchical clustering tree
)

dev.off()

write.table(file="~/reciprocal_t/analysis/GO_enrich/snp_all_GO_RESULTS", sep="\t", quote=FALSE, results)


# F1

goDatabase="go.obo"
goDivision="BP"
source("gomwu.functions.R")

input="snp_f1.txt"
goAnnotations="snp_F1_GOterms.corrected.out"

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
    perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
    largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
    smallest=5,   # a GO category should contain at least this many genes to be considered
    clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
   Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
)

pdf("~/reciprocal_t/figures/snp_F1_GO.pdf",width = 10, height = 8)

results=gomwuPlot(input,goAnnotations,goDivision,
    absValue=cutoff_snpAll,
    level1=0.05, # FDR threshold for plotting. Specify level1=1 to plot all GO categories 
    level2=0.05, # FDR cutoff to print in regular (not italic) font.
    level3=0.01,
    txtsize=1.2,    
    treeHeight=0.5, # height of the hierarchical clustering tree
)

dev.off()

write.table(file="~/reciprocal_t/analysis/GO_enrich/snp_F1_GO_RESULTS", sep="\t", quote=FALSE, results)



#######
#
# F2
#
#######

input="snp_f2.txt"
goAnnotations="snp_F2_GOterms.corrected.out"

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
    perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
    largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
    smallest=5,   # a GO category should contain at least this many genes to be considered
    clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
   Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
)

png("~/reciprocal_t/figures/snp_F2_GO.png",width = 7, height = 14, res=300, units="in")

results=gomwuPlot(input,goAnnotations,goDivision,

    absValue=cutoff_snpF2,
    level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories 
    level2=0.05, # FDR cutoff to print in regular (not italic) font.
    level3=0.01,
    txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
    treeHeight=0.5, # height of the hierarchical clustering tree
)

dev.off()

write.table(file="~/reciprocal_t/analysis/GO_enrich/snp_F2_GO_RESULTS", sep="\t", quote=FALSE, results)

#######
#
# F3
#
#######

input="snp_f3.txt"
goAnnotations="snp_F3_GOterms.corrected.out"

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
    perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
    largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
    smallest=5,   # a GO category should contain at least this many genes to be considered
    clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
   Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
)

png("~/reciprocal_t/figures/snp_F3_GO.png",width = 7, height = 10, res=300, units="in")

results=gomwuPlot(input,goAnnotations,goDivision,
    absValue=cutoff_snpF3,
    level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories 
    level2=0.05, # FDR cutoff to print in regular (not italic) font.
    level3=0.01,
    txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
    treeHeight=0.5, # height of the hierarchical clustering tree
)

dev.off()

write.table(file="~/reciprocal_t/analysis/GO_enrich/snp_F3_GO_RESULTS", sep="\t", quote=FALSE, results)

###########################
# dge
###########################


setwd("~/reciprocal_t/analysis/GO_enrich")
dgeF1 <- read.csv("~/reciprocal_t/analysis/GO_enrich/dge_f1.txt", header=FALSE, stringsAsFactors=FALSE)
dgeF2 <- read.csv("~/reciprocal_t/analysis/GO_enrich/dge_f2.txt", header=FALSE, stringsAsFactors=FALSE)
dgeF3 <- read.csv("~/reciprocal_t/analysis/GO_enrich/dge_f3.txt", header=FALSE, stringsAsFactors=FALSE)
cutoff_dgeF1 <- quantile(dgeF1$V2, 0.9, na.rm=T)
cutoff_dgeF2 <- quantile(dgeF2$V2, 0.9, na.rm=T)
cutoff_dgeF3 <- quantile(dgeF3$V2, 0.9, na.rm=T)

# summarize across all gens:
genes <- unique(c(dgeF1$V1,dgeF2$V1,dgeF3$V1))
dge_all <- data.frame(genes = genes, pval = rep(NA, length(genes)))

for(i in 1:length(genes)){
    tmpf1 <- dgeF1[which(dgeF1$V1 == genes[i]),]
    tmpf2 <- dgeF2[which(dgeF2$V1 == genes[i]),]
    tmpf3 <- dgeF3[which(dgeF3$V1 == genes[i]),]
    dge_all$pval[i] <- mean(c(tmpf1$V2,tmpf2$V2,tmpf3$V2), na.rm=T)
}

write.table(data.frame(
    x=dge_all$gene, 
    y=(dge_all$pval)), 
    file="~/reciprocal_t/analysis/GO_enrich/dge_all.txt",
    sep=",", quote=FALSE, col.names=F, row.names=F)

dgeAll <- read.csv("~/reciprocal_t/analysis/GO_enrich/dge_all.txt", header=FALSE, stringsAsFactors=FALSE)
cutoff_dgeAll <- quantile(dgeAll$V2, 0.9, na.rm=T)


goDatabase="go.obo"
goDivision="BP"
source("gomwu.functions.R")

cat ~/reciprocal_t/analysis/GO_enrich/dge_F1_GOterms.corrected.out \
    ~/reciprocal_t/analysis/GO_enrich/dge_F2_GOterms.corrected.out \
    ~/reciprocal_t/analysis/GO_enrich/dge_F3_GOterms.corrected.out |\
    sort| uniq > ~/reciprocal_t/analysis/GO_enrich/dge_ALL_GOterms.corrected.out

input="dge_all.txt"
goAnnotations="dge_ALL_GOterms.corrected.out"

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
    perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
    largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
    smallest=5,   # a GO category should contain at least this many genes to be considered
    clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
   Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
)

png("~/reciprocal_t/figures/dge_All_GO.png",width = 7, height = 14, res=300, units="in")

results=gomwuPlot(input,goAnnotations,goDivision,
    absValue=cutoff_dgeAll,
    level1=0.05, # FDR threshold for plotting. Specify level1=1 to plot all GO categories 
    level2=0.05, # FDR cutoff to print in regular (not italic) font.
    level3=0.01,
    txtsize=1.2,    
    treeHeight=0.5, # height of the hierarchical clustering tree
)

dev.off()

write.table(file="~/reciprocal_t/analysis/GO_enrich/dge_All_GO_RESULTS", sep="\t", quote=FALSE, results)


# F1


goDatabase="go.obo"
goDivision="BP"
source("gomwu.functions.R")

input="dge_f1.txt"
goAnnotations="dge_F1_GOterms.corrected.out"

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
    perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
    largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
    smallest=5,   # a GO category should contain at least this many genes to be considered
    clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
   Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
)

pdf("~/reciprocal_t/figures/dge_F1_GO.pdf",width = 10, height = 8)

results=gomwuPlot(input,goAnnotations,goDivision,
    absValue=cutoff_dgeF1,
    level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories 
    level2=0.05, # FDR cutoff to print in regular (not italic) font.
    level3=0.01,
    txtsize=1.2,    
    treeHeight=0.5, # height of the hierarchical clustering tree
)

dev.off()

write.table(file="~/reciprocal_t/analysis/GO_enrich/dge_F1_GO_RESULTS", sep="\t", quote=FALSE, results)

#######
#
# F2
#
#######

input="dge_f2.txt"
goAnnotations="dge_F2_GOterms.corrected.out"

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
    perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
    largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
    smallest=5,   # a GO category should contain at least this many genes to be considered
    clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
   Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
)

png("~/reciprocal_t/figures/dge_F2_GO.png",width = 7, height = 7, res=300, units="in")

results=gomwuPlot(input,goAnnotations,goDivision,

    absValue=cutoff_dgeF2,
    level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories 
    level2=0.05, # FDR cutoff to print in regular (not italic) font.
    level3=0.01,
    txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
    treeHeight=0.5, # height of the hierarchical clustering tree
)

dev.off()

write.table(file="~/reciprocal_t/analysis/GO_enrich/dge_F2_GO_RESULTS", sep="\t", quote=FALSE, results)

#######
#
# F3
#
#######

input="dge_f3.txt"
goAnnotations="dge_F3_GOterms.corrected.out"

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
    perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
    largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
    smallest=5,   # a GO category should contain at least this many genes to be considered
    clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
   Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
)

png("~/reciprocal_t/figures/dge_F3_GO.png",width = 7, height = 4, res=300, units="in")

results=gomwuPlot(input,goAnnotations,goDivision,
    absValue=cutoff_dgeF3,
    level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories 
    level2=0.05, # FDR cutoff to print in regular (not italic) font.
    level3=0.01,
    txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
    treeHeight=0.5, # height of the hierarchical clustering tree
)

dev.off()

write.table(file="~/reciprocal_t/analysis/GO_enrich/dge_F3_GO_RESULTS", sep="\t", quote=FALSE, results)

##########################################
####
#### pull out summary tables, format.
####
##########################################
dat <- read.csv("~/reciprocal_t/analysis/GO_enrich/MWU_BP_snp_all.txt", header=T, sep=" ")
length(which(dat$p.adj < 0.05))

filelist = list.files(pattern = "MWU_BP*")
sp_nm <- gsub("MWU_BP_","",filelist)
sp_nm <- gsub(".txt","",sp_nm)
#assuming tab separated values with a header
dfl = lapply(filelist, function(x)read.csv(x, header=T, sep=" ")) 
names(dfl) <- sp_nm


for(i in 1:length(sp_nm)){
    tmp.out <- dfl[[sp_nm[i]]][which(dfl[[sp_nm[i]]]$p.adj < 0.05),]
    write.table(file=paste("~/reciprocal_t/analysis/GO_enrich/a_", sp_nm[i], "_GO_RESULTS.txt", sep=""), 
        sep="\t", quote=FALSE, row.names=F,
        cbind(tmp.out[,c(5,6)], tmp.out[,c(2,7)]))
}

