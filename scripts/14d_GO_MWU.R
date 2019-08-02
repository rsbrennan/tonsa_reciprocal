#### determine "good genes"

snpF1 <- read.csv("~/reciprocal_t/analysis/GO_enrich/snp_meanContrib_F1.txt", header=TRUE, stringsAsFactors=FALSE)
snpF2 <- read.csv("~/reciprocal_t/analysis/GO_enrich/snp_meanContrib_F2.txt", header=TRUE, stringsAsFactors=FALSE)
snpF3 <- read.csv("~/reciprocal_t/analysis/GO_enrich/snp_meanContrib_F3.txt", header=TRUE, stringsAsFactors=FALSE)
cutoff_snpF1 <- quantile(snpF1$mean_LD1, 0.9)
cutoff_snpF2 <- quantile(snpF2$mean_LD1, 0.9)
cutoff_snpF3 <- quantile(snpF3$mean_LD1, 0.9)


goDatabase="go.obo"
goDivision="BP"
source("gomwu.functions.R")

input="snp_meanContrib_F1.txt"
goAnnotations="snp_F1_GOterms.corrected.out"

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
    perlPath="perl",
    largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
    smallest=5,   # a GO category should contain at least this many genes to be considered
    clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
   Alternative="g"
)

png("~/reciprocal_t/figures/snp_F1_GO.png",width = 7, height = 6, res=300, units="in")


results=gomwuPlot(input,goAnnotations,goDivision,
    absValue=cutoff_snpF1,
    level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories
    level2=0.05, # FDR cutoff to print in regular (not italic) font.
    level3=0.01,
    txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
    treeHeight=0.5, # height of the hierarchical clustering tree
)

dev.off()

write.table(file="~/reciprocal_t/analysis/GO_enrich/snp_F1_GO_RESULTS", sep="\t", quote=FALSE, results)

#######
#
# F2
#
#######

input="snp_meanContrib_F2.txt" # loadings scores from dapc
goAnnotations="snp_F2_GOterms.corrected.out" #snps

gomwuStats(input, goDatabase, goAnnotations, goDivision,
    perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
    largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
    smallest=5,   # a GO category should contain at least this many genes to be considered
    clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
   Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
)

png("~/reciprocal_t/figures/snp_F2_GO.png",width = 7, height = 4, res=300, units="in")

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

input="snp_meanContrib_F3.txt" # loadings scores from dapc
goAnnotations="snp_F3_GOterms.corrected.out" #snps

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
    perlPath="perl",
    largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
    smallest=5,   # a GO category should contain at least this many genes to be considered
    clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
   Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
)

png("~/reciprocal_t/figures/snp_F3_GO.png",width = 7, height = 4, res=300, units="in")

results=gomwuPlot(input,goAnnotations,goDivision,
    absValue=cutoff_snpF2,
    level1=0.1, # FDR threshold for plotting.
    level2=0.05, # FDR cutoff to print in regular (not italic) font.
    level3=0.01,
    txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
    treeHeight=0.5, # height of the hierarchical clustering tree
)

dev.off()

write.table(file="~/reciprocal_t/analysis/GO_enrich/snp_F3_GO_RESULTS", sep="\t", quote=FALSE, results)

########################################
########
######## DGE enrichment
########
########################################

dgeF1 <- read.csv("~/reciprocal_t/analysis/GO_enrich/dge_contrib_F1.go.in", header=TRUE, stringsAsFactors=FALSE)
dgeF2 <- read.csv("~/reciprocal_t/analysis/GO_enrich/dge_contrib_F2.go.in", header=TRUE, stringsAsFactors=FALSE)
dgeF3 <- read.csv("~/reciprocal_t/analysis/GO_enrich/dge_contrib_F3.go.in", header=TRUE, stringsAsFactors=FALSE)
cutoff_dgeF1 <- quantile(dgeF1$LD1, 0.9)
cutoff_dgeF2 <- quantile(dgeF2$LD1, 0.9)
cutoff_dgeF3 <- quantile(dgeF3$LD1, 0.9)


goDatabase="go.obo"
goDivision="BP"
source("gomwu.functions.R")

input="dge_contrib_F1.go.in"
goAnnotations="dge_F1_GOterms.corrected.out"

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
    perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
    largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
    smallest=5,   # a GO category should contain at least this many genes to be considered
    clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
   Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
)

png("~/reciprocal_t/figures/dge_F1_GO.png",width = 7, height = 10, res=300, units="in")

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

input="dge_contrib_F2.go.in"
goAnnotations="dge_F2_GOterms.corrected.out"

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
    perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
    largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
    smallest=5,   # a GO category should contain at least this many genes to be considered
    clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
   Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
)

png("~/reciprocal_t/figures/dge_F2_GO.png",width = 7, height = 10, res=300, units="in")

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

input="dge_contrib_F3.go.in"
goAnnotations="dge_F3_GOterms.corrected.out"

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
    perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
    largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
    smallest=5,   # a GO category should contain at least this many genes to be considered
    clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
   Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
)

png("~/reciprocal_t/figures/dge_F3_GO.png",width = 7, height = 6, res=300, units="in")

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

