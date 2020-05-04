#### determine "good genes"

########################################
########
######## DGE enrichment
########
########################################
setwd("~/reciprocal_t/analysis/GO_enrich")
dgeF1 <- read.csv("~/reciprocal_t/analysis/GO_enrich/dge_f1.txt", header=FALSE, stringsAsFactors=FALSE)
dgeF2 <- read.csv("~/reciprocal_t/analysis/GO_enrich/dge_f2.txt", header=FALSE, stringsAsFactors=FALSE)
dgeF3 <- read.csv("~/reciprocal_t/analysis/GO_enrich/dge_f3.txt", header=FALSE, stringsAsFactors=FALSE)
cutoff_dgeF1 <- quantile(dgeF1$V2, 0.9, na.rm=T)
cutoff_dgeF2 <- quantile(dgeF2$V2, 0.9, na.rm=T)
cutoff_dgeF3 <- quantile(dgeF3$V2, 0.9, na.rm=T)


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

png("~/reciprocal_t/figures/dge_F1_GO.png",width = 7, height = 14, res=300, units="in")

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
