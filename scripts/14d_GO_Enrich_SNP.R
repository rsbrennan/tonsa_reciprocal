#go enrichment for cmh test

df <- read.table("~/reciprocal_t/analysis/AF_change.txt", header=TRUE)

sigout <- df[which(df$sig_all == TRUE),]

sigout$SNP <- as.character(sigout$SNP)
df$SNP <- as.character(df$SNP)
sigout$chr <- sapply(strsplit(sigout$SNP, "\\:"), "[", 1)
df$chr <- sapply(strsplit(df$SNP, "\\:"), "[", 1)

length(unique(df$chr))

###### TOP GO #############

library(topGO)

#need to read in: cmh results, go mapping
cat ~/reciprocal_t/analysis/GO_enrich/snp_F1_GOterms.out | sed 's/\"//g' | grep -v unknown | sed 's/;/,/g' > ~/reciprocal_t/analysis/GO_enrich/snp_GOterms.nomiss.txt

sig <- read.table("~/reciprocal_t/analysis/AF_change.txt", header=TRUE)
dat <- read.csv("~/reciprocal_t/analysis/GO_enrich/snp_GOterms.nomiss.txt", header=F, sep="\t" )

sig$chr <- sapply(strsplit(as.character(sig$SNP), "\\:"), "[", 1)

#first, need a list of unique genes.
genes <- unique(sig$chr)

# also need a file that indicates the go terms for each gene, where go terms are separated by commas
ofInterest <- c()

for(i in 1:length(genes)){
    a <- sig[which(sig$chr == genes[i]),]
    if(sum(a$sig_all) > 0){
        ofInterest[i] <- TRUE
    } else{ofInterest[i] <- FALSE}
}

geneID2GO <- readMappings(file = "~/reciprocal_t/analysis/GO_enrich/snp_GOterms.nomiss.txt")
# set gene background
geneUniverse <- names(geneID2GO)

genesOfInterest <- names(geneID2GO[ofInterest])

#show genes of interest in universe vector
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

out.save <- as.data.frame(matrix(ncol=8, nrow=0))
colnames(out.save) <- c("GO.ID", "Term","Annotated","Significant","Expected","classicFisher","weight", "ontology")

for(j in c("BP", "CC", "MF")){
    myGOdata <- new("topGOdata", description="My project", 
            ontology=j, allGenes=geneList,
            annot = annFUN.gene2GO, gene2GO = geneID2GO,
            nodeSize = 5 )
        resultWeight <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
        resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
        allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, 
                    orderBy = "weight", ranksOf = "weight"  ,
                    topNodes = 200)
        allRes$ontology <- j
        out.save <- rbind(out.save,allRes[which(allRes$weight < 0.05 ),])
    }

write.table(file="~/reciprocal_t/analysis/GO_enrich/cmh_GO_enrichment.txt",
             out.save, col.names=TRUE,
            row.names=FALSE, quote=FALSE,sep="\t")


