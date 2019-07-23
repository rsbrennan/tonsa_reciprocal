


library(tximportData)
library(tximport)
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(stringr) # yes
library(reshape)
library(data.table)
library(gridExtra)


# running gene expression data first here:

#locate the directory containing the files. 
dir <- "~/reciprocal_t/analysis/salmon"
list.files(dir)

# read in table with sample ids
samples <- read.table("~/reciprocal_t/analysis/sample_id.txt", header=FALSE)

# now point to quant files
all_files <- file.path(dir, samples$V1, "quant.sf")
names(all_files) <- samples$V1

# associate transcripts with gene ids
# make a data.frame called tx2gene with two columns: 1) transcript ID and 2) gene ID. col order is important
gene_tran <- read.table("/data/atonsa/Atonsa_gen_trans_agp_gff/Atonsa_transcript_to_gene", header=FALSE)
# then convert
tx2gene <- data.frame(transcript=gene_tran$V2, gene=gene_tran$V1)

###################################################################################
###################################################################################
###
### F1
###
###################################################################################
###################################################################################

#subset files to only include F1
files <- all_files[grep("F1", all_files)]

# make sure readr package is installed
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)

# making labels for sample table
f1samp <- as.data.frame(samples$V1[grep("F1", samples$V1)])
colnames(f1samp) <- c("V1")

id <- separate(data=f1samp, col=V1, sep="_", into = c("Population", "Generation", "Replicate"))
id$Treatment <- substr(id$Population, 1,2)
id$Line <- substr(id$Population, 3,4)
id$group <- paste(id$Treatment, id$Line, sep="")

sampleTable <- data.frame(
        sampleName = f1samp$V1,
        Line = id$Line,
        Treatment = id$Treatment,
        Replicate = id$Replicate)

# assign row names to sample table
rownames(sampleTable) <- colnames(txi$counts)

# convert to DEseq
dds <- DESeqDataSetFromTximport(txi, 
                colData=sampleTable, 
                design = ~ Line + Treatment + Line:Treatment)
                      #Treatment:Line + Treatment*Line*Generation)

# remove rows where count < 10 in more than 90% of samples
keep <- apply(counts(dds), 1, function(x) {ifelse(length(which(x > 10)) > 13, TRUE, FALSE)})
dds <- dds[keep,]
nrow(dds)
#[1] 23323

# then rolog transform
rld_f1<-rlog(dds,blind=TRUE)

###################################################################################
###################################################################################
###
### F2
###
###################################################################################
###################################################################################

#subset files to only include F2
files <- all_files[grep("F2", all_files)]

# make sure readr package is installed
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)

# making labels for sample table
f2samp <- as.data.frame(samples$V1[grep("F2", samples$V1)])
colnames(f2samp) <- c("V1")

id <- separate(data=f2samp, col=V1, sep="_", into = c("Population", "Generation", "Replicate"))
id$Treatment <- substr(id$Population, 1,2)
id$Line <- substr(id$Population, 3,4)
id$group <- paste(id$Treatment, id$Line, sep="")

sampleTable <- data.frame(
        sampleName = f2samp$V1,
        Line = id$Line,
        Treatment = id$Treatment,
        Replicate = id$Replicate)

# assign row names to sample table
rownames(sampleTable) <- colnames(txi$counts)

# convert to DEseq
dds <- DESeqDataSetFromTximport(txi,
                colData=sampleTable,
                design = ~ Line + Treatment + Line:Treatment)
                      #Treatment:Line + Treatment*Line*Generation)

# remove rows where count < 10 in more than 90% of samples
keep <- apply(counts(dds), 1, function(x) {ifelse(length(which(x > 10)) > 13, TRUE, FALSE)})
dds <- dds[keep,]
nrow(dds)
#[1] 24881

# then rolog transform
rld_f2<-rlog(dds,blind=TRUE) #use blind=TRUE to not account for experimental design


###################################################################################
###################################################################################
###
### F3
###
###################################################################################
###################################################################################

#subset files to only include F1
files <- all_files[grep("F3", all_files)]

# make sure readr package is installed
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)

# making labels for sample table
f3samp <- as.data.frame(samples$V1[grep("F3", samples$V1)])
colnames(f3samp) <- c("V1")

id <- separate(data=f3samp, col=V1, sep="_", into = c("Population", "Generation", "Replicate"))
id$Treatment <- substr(id$Population, 1,2)
id$Line <- substr(id$Population, 3,4)
id$group <- paste(id$Treatment, id$Line, sep="")

sampleTable <- data.frame(
        sampleName = f3samp$V1,
        Line = id$Line,
        Treatment = id$Treatment,
        Replicate = id$Replicate)

# assign row names to sample table
rownames(sampleTable) <- colnames(txi$counts)

# convert to DEseq
dds <- DESeqDataSetFromTximport(txi, 
                colData=sampleTable, 
                design = ~ Line + Treatment + Line:Treatment)

# remove rows where count < 10 in more than 90% of samples
keep <- apply(counts(dds), 1, function(x) {ifelse(length(which(x > 10)) > 13, TRUE, FALSE)})
dds <- dds[keep,]
nrow(dds)
#[1] 24131

# then rlog transform
rld_f3<-rlog(dds,blind=TRUE) #use blind=TRUE to not account for experimental design


#######################################################
#######################################################
##
## PCA
##
#######################################################
#######################################################

data <- plotPCA(rld_f1, intgroup=c("Treatment","Line"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

data$Treatment <- gsub("AA", "AM", data$Treatment)
data$Treatment <- gsub("HH", "GH", data$Treatment)
data$Line <- gsub("AA", "AM", data$Line)
data$Line <- gsub("HH", "GH", data$Line)

a <- ggplot(data, aes(PC1, PC2, fill=Treatment, shape=Line, size=Line)) +
        geom_point() +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,24))+
        #scale_color_manual(values=c('brown3','blue3', 'darkgoldenrod3', "darkorchid3"))+
        scale_fill_manual(values=c('cornflowerblue', 'brown2'))+
        #theme(legend.position = c(0.88,0.17))+
        theme(legend.text=element_text(size=8))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
        ggtitle(" ")+
        guides(fill=guide_legend(override.aes=list(shape=21, size=7, fill=c('cornflowerblue', 'brown2')),order = 2),
                shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)),order = 1),
                size=FALSE)+
        scale_size_manual(values=c(7,5))+ theme(plot.title = element_text(hjust = 0.5))

data <- plotPCA(rld_f2, intgroup=c("Treatment","Line"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
data$Treatment <- gsub("AA", "AM", data$Treatment)
data$Treatment <- gsub("HH", "GH", data$Treatment)
data$Line <- gsub("AA", "AM", data$Line)
data$Line <- gsub("HH", "GH", data$Line)

b <- ggplot(data, aes(PC1, PC2, fill=Treatment, shape=Line, size=Line)) +
        geom_point() +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,24))+
        #scale_color_manual(values=c('brown3','blue3', 'darkgoldenrod3', "darkorchid3"))+
        scale_fill_manual(values=c('cornflowerblue', 'brown2'))+
        #theme(legend.position = c(0.88,0.17))+
        theme(legend.text=element_text(size=8))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
       ggtitle(" ")+
        guides(fill=guide_legend(override.aes=list(shape=21, size=7, fill=c('cornflowerblue', 'brown2')),order = 2),
                shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)),order = 1),
                size=FALSE)+
        scale_size_manual(values=c(7,5))+ theme(plot.title = element_text(hjust = 0.5))

data <- plotPCA(rld_f3, intgroup=c("Treatment","Line"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
data$Treatment <- gsub("AA", "AM", data$Treatment)
data$Treatment <- gsub("HH", "GH", data$Treatment)
data$Line <- gsub("AA", "AM", data$Line)
data$Line <- gsub("HH", "GH", data$Line)

c <- ggplot(data, aes(PC1, PC2, fill=Treatment, shape=Line, size=Line)) +
        geom_point() +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,24))+
        #scale_color_manual(values=c('brown3','blue3', 'darkgoldenrod3', "darkorchid3"))+
        scale_fill_manual(values=c('cornflowerblue', 'brown2'))+
        #theme(legend.position = c(0.88,0.17))+
        theme(legend.text=element_text(size=8))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
        ggtitle(" ")+
        guides(fill=guide_legend(override.aes=list(shape=21, size=7, fill=c('cornflowerblue', 'brown2')),order = 2),
                shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)),order = 1),
                size=FALSE)+
        scale_size_manual(values=c(7,5))+ theme(plot.title = element_text(hjust = 0.5))


png("~/reciprocal_t/figures/pca_RNA.png", res=300, height=4, width=11, units="in")

ggarrange(a, b,c nrow = 1, ncol=3, common.legend=TRUE)

dev.off()


##### snps


af <- read.table("~/reciprocal_t/analysis/filtered_allele_freqs.txt", header=TRUE)
dat3 <- read.table("~/reciprocal_t/analysis/filtered_variants.txt", header=TRUE)

pops <- c(
        "AAAA_F1_REP1", "AAAA_F1_REP2", "AAAA_F1_REP3", "AAAA_F1_REP4",
        "AAAA_F2_REP1", "AAAA_F2_REP2", "AAAA_F2_REP3", "AAAA_F2_REP4",
        "AAAA_F3_REP1", "AAAA_F3_REP2", "AAAA_F3_REP3", "AAAA_F3_REP4",
        "AAHH_F1_REP1", "AAHH_F1_REP2", "AAHH_F1_REP3", "AAHH_F1_REP4",
        "AAHH_F2_REP1", "AAHH_F2_REP2", "AAHH_F2_REP3", "AAHH_F2_REP4",
        "AAHH_F3_REP1", "AAHH_F3_REP2", "AAHH_F3_REP3", "AAHH_F3_REP4",
        "HHAA_F1_REP1", "HHAA_F1_REP2", "HHAA_F1_REP3", "HHAA_F1_REP4",
        "HHAA_F2_REP1", "HHAA_F2_REP2", "HHAA_F2_REP3", "HHAA_F2_REP4",
        "HHAA_F3_REP1", "HHAA_F3_REP2", "HHAA_F3_REP3", "HHAA_F3_REP4",
        "HHHH_F1_REP1", "HHHH_F1_REP2", "HHHH_F1_REP3", "HHHH_F1_REP4",
        "HHHH_F2_REP1", "HHHH_F2_REP2", "HHHH_F2_REP3", "HHHH_F2_REP4",
        "HHHH_F3_REP1", "HHHH_F3_REP2", "HHHH_F3_REP3", "HHHH_F3_REP4")


freqs <- t(af[,2:ncol(af)])
colnames(freqs) <- c(paste(dat3$Chrom, dat3$Position, sep=":"))

####
##
## plot pca
##
####

# f1 alone

freqs_f1 <- freqs[grep("F1", row.names(freqs)),]
freqs_f1 <- freqs_f1[,which(apply(freqs_f1, 2, var)>0)]
pcaResult <- prcomp(freqs_f1, scale=TRUE)

percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

f1 <- pops[grep("F1", pops)]

data <- data.frame(id=f1, Line=substr(f1, 3,4), 
    Treatment=substr(f1, 1,2),
        gen=substr(f1, 6,7),
        PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,2])

data$Treatment <- gsub("AA", "AM", data$Treatment)
data$Treatment <- gsub("HH", "GH", data$Treatment)
data$Line <- gsub("AA", "AM", data$Line)
data$Line <- gsub("HH", "GH", data$Line)

d <- ggplot(data, aes(PC1, PC2, fill=Treatment, shape=Line, size=Line)) +
        geom_point() +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,24))+
        #scale_color_manual(values=c('brown3','blue3', 'darkgoldenrod3', "darkorchid3"))+
        scale_fill_manual(values=c('cornflowerblue', 'brown2'))+
        #theme(legend.position = c(0.88,0.17))+
        theme(legend.text=element_text(size=8))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
        ggtitle("F1")+
        guides(fill=guide_legend(override.aes=list(shape=21, size=7, fill=c('cornflowerblue', 'brown2')),order = 2),
                shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)),order = 1),
                size=FALSE)+
        scale_size_manual(values=c(7,5)) + theme(plot.title = element_text(hjust = 0.5))


# f2

freqs_f2 <- freqs[grep("F2", row.names(freqs)),]
freqs_f2 <- freqs_f2[,which(apply(freqs_f2, 2, var)>0)]
pcaResult <- prcomp(freqs_f2, scale=TRUE)

percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

f2 <- pops[grep("F2", pops)]

data <- data.frame(id=f2, Line=substr(f2, 3,4),
    Treatment=substr(f2, 1,2),
        gen=substr(f2, 6,7),
        PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,2])
data$Treatment <- gsub("AA", "AM", data$Treatment)
data$Treatment <- gsub("HH", "GH", data$Treatment)
data$Line <- gsub("AA", "AM", data$Line)
data$Line <- gsub("HH", "GH", data$Line)

e <- ggplot(data, aes(PC1, PC2, fill=Treatment, shape=Line, size=Line)) +
        geom_point() +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,24))+
        #scale_color_manual(values=c('brown3','blue3', 'darkgoldenrod3', "darkorchid3"))+
        scale_fill_manual(values=c('cornflowerblue', 'brown2'))+
        #theme(legend.position = c(0.88,0.17))+
        theme(legend.text=element_text(size=8))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
        ggtitle("F2")+
        guides(fill=guide_legend(override.aes=list(shape=21, size=7, fill=c('cornflowerblue', 'brown2')),order = 2),
                shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)),order = 1),
                size=FALSE)+
        scale_size_manual(values=c(7,5)) +theme(plot.title = element_text(hjust = 0.5))

# f3

freqs_f3 <- freqs[grep("F3", row.names(freqs)),]
freqs_f3 <- freqs_f3[,which(apply(freqs_f3, 2, var)>0)]
pcaResult <- prcomp(freqs_f3, scale=TRUE)

percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

f3 <- pops[grep("F3", pops)]

data <- data.frame(id=f3, Line=substr(f3, 3,4), 
    Treatment=substr(f3, 1,2),
        gen=substr(f3, 6,7),
        PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,2])
data$Treatment <- gsub("AA", "AM", data$Treatment)
data$Treatment <- gsub("HH", "GH", data$Treatment)
data$Line <- gsub("AA", "AM", data$Line)
data$Line <- gsub("HH", "GH", data$Line)

f <- ggplot(data, aes(PC1, PC2, fill=Treatment, shape=Line, size=Line)) +
        geom_point() +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,24))+
        #scale_color_manual(values=c('brown3','blue3', 'darkgoldenrod3', "darkorchid3"))+
        scale_fill_manual(values=c('cornflowerblue', 'brown2'))+
        #theme(legend.position = c(0.88,0.17))+
        theme(legend.text=element_text(size=8))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
        ggtitle("F3")+
        guides(fill=guide_legend(override.aes=list(shape=21, size=7, fill=c('cornflowerblue', 'brown2')),order = 2),
                shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)),order = 1),
                size=FALSE)+
        scale_size_manual(values=c(7,5)) + theme(plot.title = element_text(hjust = 0.5))


png("~/reciprocal_t/figures/pca_SNP.png", res=300, height=4, width=11, units="in")

ggarrange(b, c, d, nrow = 1, ncol=3, common.legend=TRUE)

dev.off()

pdf("~/reciprocal_t/figures/fig1_pca.pdf",height=8, width=11)

ggarrange(d,e,f,a,b,c, nrow = 2, ncol=3, common.legend=TRUE, labels="AUTO", legend="bottom")

dev.off()
