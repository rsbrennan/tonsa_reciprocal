
# script for correlation plots for DGE
## this plots vertical, not horizontal

library(tximport)
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(scales)

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
gene_tran <- read.table("/data/copepods/tonsa_transcriptome/Atonsa_gen_trans_agp_gff/Atonsa_transcript_to_gene", header=FALSE)
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

dat_f1=as.data.frame(assay(rld_f1))
colnames(dat_f1)<-colnames(dds)

write.table(dat_f1, file="~/reciprocal_t/analysis/F1_rlog_counts.txt", sep="\t", quote=FALSE)

##########################
#### run DEseq
##########################

# change levels. so I know what the reference level is. 
dds$Treatment <- relevel(dds$Treatment, ref = "AA")
dds$Line <- relevel(dds$Line, ref = "AA")
dds$group <- factor(paste0(dds$Treatment, dds$Line))

design(dds) <- ~ group

dds$group <- relevel(dds$group, ref = "AAAA")

dds_f1 <- DESeq(dds, parallel=T)
res_f1 <- results(dds_f1, alpha=0.05)

dds <- estimateSizeFactors(dds)

write.table(counts(dds, normalized=TRUE), file="~/reciprocal_t/analysis/F1_norm_counts.txt", sep="\t", quote=FALSE)

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

dat_f2=as.data.frame(assay(rld_f2))
colnames(dat_f2)<-colnames(dds)

write.table(dat_f2, file="~/reciprocal_t/analysis/F2_rlog_counts.txt", sep="\t", quote=FALSE) 

##########################
#### run DEseq
##########################

# change levels. so I know what the reference level is.
dds$Treatment <- relevel(dds$Treatment, ref = "AA")
dds$Line <- relevel(dds$Line, ref = "AA")
dds$group <- factor(paste0(dds$Treatment, dds$Line))

design(dds) <- ~ group
dds$group <- relevel(dds$group, ref = "AAAA")

dds_f2 <- DESeq(dds, parallel=T)
res_f2 <- results(dds_f2, alpha=0.05)

dds <- estimateSizeFactors(dds)

write.table(counts(dds, normalized=TRUE), file="~/reciprocal_t/analysis/F2_norm_counts.txt", sep="\t", quote=FALSE)


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

dat_f3=as.data.frame(assay(rld_f3)) 
colnames(dat_f3)<-colnames(dds)

write.table(dat_f3, file="~/reciprocal_t/analysis/F3_rlog_counts.txt", sep="\t", quote=FALSE) 

##########################
#### run DEseq
##########################

# change levels. so I know what the reference level is. 
dds$Treatment <- relevel(dds$Treatment, ref = "AA")
dds$Line <- relevel(dds$Line, ref = "AA")
dds$group <- factor(paste0(dds$Treatment, dds$Line))

design(dds) <- ~ group
dds$group <- relevel(dds$group, ref = "AAAA")

dds_f3 <- DESeq(dds, parallel=T)
res_f3 <- results(dds_f3, alpha=0.05)

dds <- estimateSizeFactors(dds)

write.table(counts(dds, normalized=TRUE), file="~/reciprocal_t/analysis/F3_norm_counts.txt", sep="\t", quote=FALSE)


########################################################################
####################################
######
######  adaptive change in expression
######
####################################
########################################################################

# pull out adaptively diverged transcripts
#######
### HHHH vs AAAA
#######

# AA plasticity

res_HHHHvsAAAA_f1 <- results(dds_f1, contrast=c("group", "HHHH", "AAAA"))
# summary(res_HHHHvsAAAA_f1, alpha=0.05)
# sum(res_HHHHvsAAAA_f1$padj < 0.05, na.rm=TRUE)

# parse down to the ones that show adaptive differentiation
sig <- as.data.frame(res_HHHHvsAAAA_f1)[which(as.data.frame(res_HHHHvsAAAA_f1$padj) < 0.05),]

res_AAAAvsHHAA_f1 <- results(dds_f1, contrast=c("group", "HHAA", "AAAA"))
#summary(res_AAAAvsHHAA_f1, alpha=0.05)
#sum(res_AAAAvsHHAA_f1$padj < 0.05, na.rm=TRUE)

f1_AA_plasticity <- merge(as.data.frame(res_AAAAvsHHAA_f1), sig, by="row.names")

names(f1_AA_plasticity)[names(f1_AA_plasticity) == 'log2FoldChange.x'] <- 'log2FoldChange_AAAAvsHHAA'
names(f1_AA_plasticity)[names(f1_AA_plasticity) == 'log2FoldChange.y'] <- 'log2FoldChange_HHHHvsAAAA'

f1_AA_out <- (merge(as.data.frame(res_AAAAvsHHAA_f1), as.data.frame(res_HHHHvsAAAA_f1),
                    by="row.names"))

# f2

res_HHHHvsAAAA_f2 <- results(dds_f2, contrast=c("group", "HHHH", "AAAA"))
#summary(res_HHHHvsAAAA_f2, alpha=0.05)
#sum(res_HHHHvsAAAA_f2$padj < 0.05, na.rm=TRUE)

# parse down to the ones that show adaptive differentiation
sig <- as.data.frame(res_HHHHvsAAAA_f2)[which(as.data.frame(res_HHHHvsAAAA_f2$padj) < 0.05),]

res_AAAAvsHHAA_f2 <- results(dds_f2, contrast=c("group", "HHAA", "AAAA"))
#summary(res_AAAAvsHHAA_f2, alpha=0.05)
#sum(res_AAAAvsHHAA_f2$padj < 0.05, na.rm=TRUE)

f2_AA_plasticity <- merge(as.data.frame(res_AAAAvsHHAA_f2), sig, by="row.names")

names(f2_AA_plasticity)[names(f2_AA_plasticity) == 'log2FoldChange.x'] <- 'log2FoldChange_AAAAvsHHAA'
names(f2_AA_plasticity)[names(f2_AA_plasticity) == 'log2FoldChange.y'] <- 'log2FoldChange_HHHHvsAAAA'

f2_AA_out <- (merge(as.data.frame(res_AAAAvsHHAA_f2), as.data.frame(res_HHHHvsAAAA_f2),
                    by="row.names"))

# f3

res_HHHHvsAAAA_f3 <- results(dds_f3, contrast=c("group", "HHHH", "AAAA"))
#summary(res_HHHHvsAAAA_f3, alpha=0.05)
#sum(res_HHHHvsAAAA_f3$padj < 0.05, na.rm=TRUE)

# parse down to the ones that show adaptive differentiation
sig <- as.data.frame(res_HHHHvsAAAA_f3)[which(as.data.frame(res_HHHHvsAAAA_f3$padj) < 0.05),]

res_AAAAvsHHAA_f3 <- results(dds_f3, contrast=c("group", "HHAA", "AAAA"))
#summary(res_AAAAvsHHAA_f3, alpha=0.05)
#sum(res_AAAAvsHHAA_f3$padj < 0.05, na.rm=TRUE)

f3_AA_plasticity <- merge(as.data.frame(res_AAAAvsHHAA_f3), sig, by="row.names")

names(f3_AA_plasticity)[names(f3_AA_plasticity) == 'log2FoldChange.x'] <- 'log2FoldChange_AAAAvsHHAA'
names(f3_AA_plasticity)[names(f3_AA_plasticity) == 'log2FoldChange.y'] <- 'log2FoldChange_HHHHvsAAAA'

f3_AA_out <- (merge(as.data.frame(res_AAAAvsHHAA_f3), as.data.frame(res_HHHHvsAAAA_f3),
                    by="row.names"))

#######
## HH plasticity
#######

res_HHHHvsAAAA_f1 <- results(dds_f1, contrast=c("group", "AAAA", "HHHH"))

# parse down to the ones that show adaptive differentiation
sig <- as.data.frame(res_HHHHvsAAAA_f1)[which(as.data.frame(res_HHHHvsAAAA_f1$padj) < 0.05),]

res_HHHHvsAAHH_f1 <- results(dds_f1, contrast=c("group", "AAHH", "HHHH"))
#summary(res_HHHHvsAAHH_f1, alpha=0.05)
#sum(res_HHHHvsAAHH_f1$padj < 0.05, na.rm=TRUE)

f1_HH_plasticity <- merge(as.data.frame(res_HHHHvsAAHH_f1), sig, by="row.names")

names(f1_HH_plasticity)[names(f1_HH_plasticity) == 'log2FoldChange.x'] <- 'log2FoldChange_HHHHvsAAHH'
names(f1_HH_plasticity)[names(f1_HH_plasticity) == 'log2FoldChange.y'] <- 'log2FoldChange_HHHHvsAAAA'

f1_HH_out <- (merge(as.data.frame(res_HHHHvsAAHH_f1), as.data.frame(res_HHHHvsAAAA_f1),
                    by="row.names"))

# F2

res_HHHHvsAAAA_f2 <- results(dds_f2, contrast=c("group", "AAAA", "HHHH"))

# parse down to the ones that show adaptive differentiation
sig <- as.data.frame(res_HHHHvsAAAA_f2)[which(as.data.frame(res_HHHHvsAAAA_f2$padj) < 0.05),]

res_HHHHvsAAHH_f2 <- results(dds_f2, contrast=c("group", "AAHH", "HHHH"))
#summary(res_HHHHvsAAHH_f2, alpha=0.05)
#sum(res_HHHHvsAAHH_f2$padj < 0.05, na.rm=TRUE)

f2_HH_plasticity <- merge(as.data.frame(res_HHHHvsAAHH_f2), sig, by="row.names")

names(f2_HH_plasticity)[names(f2_HH_plasticity) == 'log2FoldChange.x'] <- 'log2FoldChange_HHHHvsAAHH'
names(f2_HH_plasticity)[names(f2_HH_plasticity) == 'log2FoldChange.y'] <- 'log2FoldChange_HHHHvsAAAA'

f2_HH_out <- (merge(as.data.frame(res_HHHHvsAAHH_f2), as.data.frame(res_HHHHvsAAAA_f2),
                    by="row.names"))

# F3

res_HHHHvsAAAA_f3 <- results(dds_f3, contrast=c("group", "AAAA", "HHHH"))

# parse down to the ones that show adaptive differentiation
sig <- as.data.frame(res_HHHHvsAAAA_f3)[which(as.data.frame(res_HHHHvsAAAA_f3$padj) < 0.05),]

res_HHHHvsAAHH_f3 <- results(dds_f3, contrast=c("group", "AAHH", "HHHH"))
#summary(res_HHHHvsAAHH_f3, alpha=0.05)
#sum(res_HHHHvsAAHH_f3$padj < 0.05, na.rm=TRUE)

f3_HH_plasticity <- merge(as.data.frame(res_HHHHvsAAHH_f3), sig, by="row.names")

names(f3_HH_plasticity)[names(f3_HH_plasticity) == 'log2FoldChange.x'] <- 'log2FoldChange_HHHHvsAAHH'
names(f3_HH_plasticity)[names(f3_HH_plasticity) == 'log2FoldChange.y'] <- 'log2FoldChange_HHHHvsAAAA'

f3_HH_out <- (merge(as.data.frame(res_HHHHvsAAHH_f3), as.data.frame(res_HHHHvsAAAA_f3),
                    by="row.names"))
########################################################################
####################################
######  Plotting
####################################
######################################################################## 

pdf("~/reciprocal_t/figures/DGE_AM.pdf", height = 3.75, width = 1.75)

par(mfrow = c(3, 1))
par(cex = 0.6)
par(mar = c(0, 3.7, 0, 0), oma = c(3.8, 0, 0.5, 0.9))
par(tcl = -0.25)
par(mgp = c(1.5, 0.3, 0))

plot(x= f1_AA_plasticity$log2FoldChange_HHHHvsAAAA, 
    y=f1_AA_plasticity$log2FoldChange_AAAAvsHHAA,
    ylab="", 
    xlab= "",
    col=alpha("black", 0.4), bg=alpha("gray38", 1), pch=21, 
    ylim=c(-6.5, 6.5), xlim=c(-9, 9),
    cex.lab=1, cex.main=1.4, cex=1.2, xaxt='n', yaxt='n')
axis(2, at=c(-6, -3, 0, 3, 6), cex.axis=0.8)

nonAdapt_f1_AA <- f1_AA_plasticity[which(
    (f1_AA_plasticity$log2FoldChange_AAAAvsHHAA > 0 & f1_AA_plasticity$log2FoldChange_HHHHvsAAAA < 0) | 
    (f1_AA_plasticity$log2FoldChange_AAAAvsHHAA < 0 & f1_AA_plasticity$log2FoldChange_HHHHvsAAAA > 0)),]
points(x= nonAdapt_f1_AA$log2FoldChange_HHHHvsAAAA, y=nonAdapt_f1_AA$log2FoldChange_AAAAvsHHAA,
    col=alpha("black", 0.4), bg=alpha("gray86", 1), pch=21, cex=1.2)
Adapt_f1_AA <- f1_AA_plasticity[which(
    (f1_AA_plasticity$log2FoldChange_AAAAvsHHAA < 0 & f1_AA_plasticity$log2FoldChange_HHHHvsAAAA < 0) |
    (f1_AA_plasticity$log2FoldChange_AAAAvsHHAA > 0 & f1_AA_plasticity$log2FoldChange_HHHHvsAAAA > 0)),]

abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(f1_AA_plasticity$log2FoldChange_AAAAvsHHAA ~ f1_AA_plasticity$log2FoldChange_HHHHvsAAAA), 
col="firebrick3", lwd=1.5)

#summary(lm(f1_AA_plasticity$log2FoldChange_AAAAvsHHAA ~ f1_AA_plasticity$log2FoldChange_HHHHvsAAAA))

corval <- round(cor.test(f1_AA_plasticity$log2FoldChange_HHHHvsAAAA,f1_AA_plasticity$log2FoldChange_AAAAvsHHAA, 
    method="pearson")$estimate, 2)
text(x=-5, y=4, bquote(rho == .(corval)), cex=1.2)


# AA f2

plot(x= f2_AA_plasticity$log2FoldChange_HHHHvsAAAA, y=f2_AA_plasticity$log2FoldChange_AAAAvsHHAA,
    ylab="", 
    xlab= "",
    col=alpha("black", 0.4), bg=alpha("gray38", 1), pch=21, 
    ylim=c(-6.5, 6.5), xlim=c(-9, 9),
    cex.lab=1, cex.main=1.4, cex=1.2, xaxt='n', yaxt='n')
axis(2, at=c(-6, -3, 0, 3, 6), cex.axis=0.8)

nonAdapt_f2_AA <- f2_AA_plasticity[which(
    (f2_AA_plasticity$log2FoldChange_AAAAvsHHAA > 0 & f2_AA_plasticity$log2FoldChange_HHHHvsAAAA < 0) | 
    (f2_AA_plasticity$log2FoldChange_AAAAvsHHAA < 0 & f2_AA_plasticity$log2FoldChange_HHHHvsAAAA > 0)),]
points(x= nonAdapt_f2_AA$log2FoldChange_HHHHvsAAAA, y=nonAdapt_f2_AA$log2FoldChange_AAAAvsHHAA,
    col=alpha("black", 0.4), bg=alpha("gray86", 1), pch=21, cex=1.2)
Adapt_f2_AA <- f2_AA_plasticity[which(
     (f2_AA_plasticity$log2FoldChange_AAAAvsHHAA < 0 & f2_AA_plasticity$log2FoldChange_HHHHvsAAAA < 0) |
     (f2_AA_plasticity$log2FoldChange_AAAAvsHHAA > 0 & f2_AA_plasticity$log2FoldChange_HHHHvsAAAA > 0)),]


abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(f2_AA_plasticity$log2FoldChange_AAAAvsHHAA~ f2_AA_plasticity$log2FoldChange_HHHHvsAAAA), 
    col="firebrick3", lwd=1.5)

#summary(lm(f2_AA_plasticity$log2FoldChange_AAAAvsHHAA~ f2_AA_plasticity$log2FoldChange_HHHHvsAAAA))

corval <- round(cor.test(f2_AA_plasticity$log2FoldChange_HHHHvsAAAA,f2_AA_plasticity$log2FoldChange_AAAAvsHHAA, 
    method="pearson")$estimate, 2)
text(x=-5, y=4, bquote(rho == .(corval)), cex=1.2)

 par(xpd=TRUE)
mtext(text = "Plastic change in expression\nlog2 fold change: AM_in_GH/AM",
      side = 2,#side 1 = bottom
      line = 1.7, cex=0.6)
 par(xpd=FALSE)

# AA f3

plot(x= f3_AA_plasticity$log2FoldChange_HHHHvsAAAA, y=f3_AA_plasticity$log2FoldChange_AAAAvsHHAA,
    ylab="", 
    xlab= "",
    col=alpha("black", 0.4), bg=alpha("gray38", 1), pch=21, 
    ylim=c(-6.5, 6.5), xlim=c(-9, 9),
    cex.lab=1, cex.main=1.4, cex=1.2, xaxt='n', yaxt='n')
axis(2, at=c(-6, -3, 0, 3, 6), cex.axis=0.8)
axis(1, cex.axis=0.8)

nonAdapt_f3_AA <- f3_AA_plasticity[which(
    (f3_AA_plasticity$log2FoldChange_AAAAvsHHAA > 0 & f3_AA_plasticity$log2FoldChange_HHHHvsAAAA < 0) | 
    (f3_AA_plasticity$log2FoldChange_AAAAvsHHAA < 0 & f3_AA_plasticity$log2FoldChange_HHHHvsAAAA > 0)),]
points(x= nonAdapt_f3_AA$log2FoldChange_HHHHvsAAAA, y=nonAdapt_f3_AA$log2FoldChange_AAAAvsHHAA,
    col=alpha("black", 0.4), bg=alpha("gray86", 1), pch=21, cex=1.2)
Adapt_f3_AA <- f3_AA_plasticity[which(
     (f3_AA_plasticity$log2FoldChange_AAAAvsHHAA < 0 & f3_AA_plasticity$log2FoldChange_HHHHvsAAAA < 0) |
     (f3_AA_plasticity$log2FoldChange_AAAAvsHHAA > 0 & f3_AA_plasticity$log2FoldChange_HHHHvsAAAA > 0)),]

abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(f3_AA_plasticity$log2FoldChange_AAAAvsHHAA~ f3_AA_plasticity$log2FoldChange_HHHHvsAAAA), col="firebrick3", 
    lwd=1.5)

#summary(lm(f3_AA_plasticity$log2FoldChange_AAAAvsHHAA~ f3_AA_plasticity$log2FoldChange_HHHHvsAAAA))

corval <- round(cor.test(f3_AA_plasticity$log2FoldChange_HHHHvsAAAA,f3_AA_plasticity$log2FoldChange_AAAAvsHHAA, 
    method="pearson")$estimate, 2)
text(x=-5, y=4, bquote(rho == .(corval)), cex=1.2)

 par(xpd=TRUE)
mtext(text = "Evolved change in expression:\nlog2 fold change: GH/AM",
      side = 1,#side 1 = bottom
      line = 2.3, cex=0.55)
 par(xpd=FALSE)


dev.off()



################
#### HH
################

pdf("~/reciprocal_t/figures/DGE_GH.pdf",height = 3.75, width = 1.75)

par(mfrow = c(3, 1))
par(mar = c(0, 3.7, 0, 0), oma = c(3.8, 0, 0.5, 0.9))
par(tcl = -0.25)
par(cex = 0.6)
par(mgp = c(1.5, 0.3, 0))

plot(x= f1_HH_plasticity$log2FoldChange_HHHHvsAAAA, y=f1_HH_plasticity$log2FoldChange_HHHHvsAAHH,
    ylab="",
    xlab= "",
    col=alpha("black", 0.4), bg=alpha("gray38", 1), pch=21,
    ylim=c(-6.5, 6.5), xlim=c(-9, 9),
    cex.lab=1, cex.main=1.4, cex=1.2, xaxt='n', yaxt='n')
axis(2, at=c(-6, -3, 0, 3, 6), cex.axis=0.8)

nonAdapt_f1_HH <- f1_HH_plasticity[which(
    (f1_HH_plasticity$log2FoldChange_HHHHvsAAHH > 0 & f1_HH_plasticity$log2FoldChange_HHHHvsAAAA < 0) |
    (f1_HH_plasticity$log2FoldChange_HHHHvsAAHH < 0 & f1_HH_plasticity$log2FoldChange_HHHHvsAAAA > 0)),]
points(x= nonAdapt_f1_HH$log2FoldChange_HHHHvsAAAA, y=nonAdapt_f1_HH$log2FoldChange_HHHHvsAAHH,
    col=alpha("black", 0.4), bg=alpha("gray86", 1), pch=21, cex=1.2)

Adapt_f1_HH <- f1_HH_plasticity[which(
    (f1_HH_plasticity$log2FoldChange_HHHHvsAAHH < 0 & f1_HH_plasticity$log2FoldChange_HHHHvsAAAA < 0) |
    (f1_HH_plasticity$log2FoldChange_HHHHvsAAHH > 0 & f1_HH_plasticity$log2FoldChange_HHHHvsAAAA > 0)),]

abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(f1_HH_plasticity$log2FoldChange_HHHHvsAAHH~ f1_HH_plasticity$log2FoldChange_HHHHvsAAAA), 
    col="firebrick3", lwd=1.5)

corval <- round(cor.test(f1_HH_plasticity$log2FoldChange_HHHHvsAAAA,f1_HH_plasticity$log2FoldChange_HHHHvsAAHH, 
    method="pearson")$estimate, 2)
text(x=-5, y=4, bquote(rho == .(corval)), cex=1.2)

#### HH F2

plot(x= f2_HH_plasticity$log2FoldChange_HHHHvsAAAA, y=f2_HH_plasticity$log2FoldChange_HHHHvsAAHH,
    ylab="", 
    xlab= "",
    col=alpha("black", 0.4), bg=alpha("gray38", 1), pch=21, 
    ylim=c(-6.5, 6.5), xlim=c(-9, 9),
    cex.lab=1, cex.main=1.4, cex=1.2, xaxt='n', yaxt='n')
axis(2, at=c(-6, -3, 0, 3, 6), cex.axis=0.8)

nonAdapt_f2_HH <- f2_HH_plasticity[which(
    (f2_HH_plasticity$log2FoldChange_HHHHvsAAHH > 0 & f2_HH_plasticity$log2FoldChange_HHHHvsAAAA < 0) | 
    (f2_HH_plasticity$log2FoldChange_HHHHvsAAHH < 0 & f2_HH_plasticity$log2FoldChange_HHHHvsAAAA > 0)),]
points(x= nonAdapt_f2_HH$log2FoldChange_HHHHvsAAAA, y=nonAdapt_f2_HH$log2FoldChange_HHHHvsAAHH,
    col=alpha("black", 0.4), bg=alpha("gray86", 1), pch=21, cex=1.2)
Adapt_f2_HH <- f2_HH_plasticity[which(
     (f2_HH_plasticity$log2FoldChange_HHHHvsAAHH < 0 & f2_HH_plasticity$log2FoldChange_HHHHvsAAAA < 0) |
     (f2_HH_plasticity$log2FoldChange_HHHHvsAAHH > 0 & f2_HH_plasticity$log2FoldChange_HHHHvsAAAA > 0)),]

abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(f2_HH_plasticity$log2FoldChange_HHHHvsAAHH ~ f2_HH_plasticity$log2FoldChange_HHHHvsAAAA), 
    col="firebrick3", lwd=1.5)

corval <- round(cor.test(f2_HH_plasticity$log2FoldChange_HHHHvsAAAA,f2_HH_plasticity$log2FoldChange_HHHHvsAAHH, 
    method="pearson")$estimate, 2)
text(x=-5, y=4, bquote(rho == .(corval)), cex=1.2)


 par(xpd=TRUE)
mtext(text = "Plastic change in expression\nlog2 fold change: GH_in_AM/GH",
      side = 2,#side 1 = bottom
      line = 1.7, cex=0.6)
 par(xpd=FALSE)


#### HH F3

plot(x= f3_HH_plasticity$log2FoldChange_HHHHvsAAAA, y=f3_HH_plasticity$log2FoldChange_HHHHvsAAHH,
    ylab="", 
    xlab= "",
    col=alpha("black", 0.4), bg=alpha("gray38", 1), pch=21, 
    ylim=c(-6.5, 6.5), xlim=c(-9, 9),
    cex.lab=1, cex.main=1.4, cex=1.2, xaxt='n', yaxt='n')
axis(2, at=c(-6, -3, 0, 3, 6), cex.axis=0.8)
axis(1, cex.axis=0.8)

nonAdapt_f3_HH <- f3_HH_plasticity[which(
    (f3_HH_plasticity$log2FoldChange_HHHHvsAAHH > 0 & f3_HH_plasticity$log2FoldChange_HHHHvsAAAA < 0) | 
    (f3_HH_plasticity$log2FoldChange_HHHHvsAAHH < 0 & f3_HH_plasticity$log2FoldChange_HHHHvsAAAA > 0)),]
points(x= nonAdapt_f3_HH$log2FoldChange_HHHHvsAAAA, y=nonAdapt_f3_HH$log2FoldChange_HHHHvsAAHH,
    col=alpha("black", 0.4), bg=alpha("gray86", 1), pch=21, cex=1.2)
Adapt_f3_HH <- f3_HH_plasticity[which(
     (f3_HH_plasticity$log2FoldChange_HHHHvsAAHH < 0 & f3_HH_plasticity$log2FoldChange_HHHHvsAAAA < 0) |
     (f3_HH_plasticity$log2FoldChange_HHHHvsAAHH > 0 & f3_HH_plasticity$log2FoldChange_HHHHvsAAAA > 0)),]

abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(f3_HH_plasticity$log2FoldChange_HHHHvsAAHH ~ f3_HH_plasticity$log2FoldChange_HHHHvsAAAA), 
    col="firebrick3", lwd=1.5)


corval <- round(cor.test(f3_HH_plasticity$log2FoldChange_HHHHvsAAAA,f3_HH_plasticity$log2FoldChange_HHHHvsAAHH, 
    method="pearson")$estimate, 2)
text(x=-5, y=4, bquote(rho == .(corval)), cex=1.2)


 par(xpd=TRUE)
mtext(text = "Evolved change in expression:\nlog2 fold change: AM/GH",
      side = 1,#side 1 = bottom
      line = 2.3, cex=0.55)
 par(xpd=FALSE)

dev.off()


# prop adaptive
(nrow(f1_HH_plasticity)-nrow(Adapt_f1_HH))/nrow(f1_HH_plasticity)
# [1] 0.3457014
(nrow(f1_AA_plasticity)-nrow(Adapt_f1_AA))/nrow(f1_AA_plasticity)
# [1] 0.1466063
(nrow(f2_HH_plasticity)-nrow(Adapt_f2_HH))/nrow(f2_HH_plasticity)
# [1] 0.05859375
(nrow(f2_AA_plasticity)-nrow(Adapt_f2_AA))/nrow(f2_AA_plasticity)
# [1] 0.1269531
(nrow(f3_HH_plasticity)-nrow(Adapt_f3_HH))/nrow(f3_HH_plasticity)
# [1] 0.05309735
(nrow(f3_AA_plasticity)-nrow(Adapt_f3_AA))/nrow(f3_AA_plasticity)
# [1] 0.1017699

prop.test(c(nrow(Adapt_f1_HH),nrow(Adapt_f2_HH),nrow(Adapt_f3_HH)),
          c(nrow(f1_HH_plasticity),nrow(f2_HH_plasticity),nrow(f3_HH_plasticity)),correct=FALSE)
# X-squared = 259.81, df = 2, p-value < 2.2e-16
prop.test(c(nrow(Adapt_f2_HH),nrow(Adapt_f3_HH)),c(nrow(f2_HH_plasticity),nrow(f3_HH_plasticity)),correct=FALSE)
# post hoc shows that F1 is diff from F2 and F3

prop.test(c(nrow(Adapt_f1_AA),nrow(Adapt_f2_AA),nrow(Adapt_f3_AA)),
          c(nrow(f1_AA_plasticity),nrow(f2_AA_plasticity),nrow(f3_AA_plasticity)),correct=FALSE)
# X-squared = 5.7784, df = 2, p-value = 0.05562

prop.test(c(nrow(Adapt_f1_AA),nrow(Adapt_f2_AA)),
          c(nrow(f1_AA_plasticity),nrow(f2_AA_plasticity)),correct=FALSE)
# X-squared = 1.1199, df = 1, p-value = 0.2899
prop.test(c(nrow(Adapt_f3_AA),nrow(Adapt_f1_AA)),
          c(nrow(f3_AA_plasticity),nrow(f1_AA_plasticity)),correct=FALSE)
# X-squared = 5.5715, df = 1, p-value = 0.01825

# note that the prop.test is just a chisq test. the same as folling:
#chisq.test(matrix(c(nrow(Adapt_f3_AA),nrow(Adapt_f1_AA),
#         nrow(nonAdapt_f3_AA),nrow(nonAdapt_f1_AA)), nrow=2, ncol=2),correct=FALSE)

# and that an exact test gives more or less the same results
#fisher.test(matrix(c(nrow(Adapt_f3_AA),nrow(Adapt_f1_AA),
#         nrow(nonAdapt_f3_AA),nrow(nonAdapt_f1_AA)), nrow=2, ncol=2))





####
#
# save output from results
#
####

write.table(data.frame(x=row.names(as.data.frame(res_HHHHvsAAAA_f1)), y=-log10(as.data.frame(res_HHHHvsAAAA_f1)[,5])), 
    file="~/reciprocal_t/analysis/GO_enrich/dge_f1.txt",
    sep=",", quote=FALSE, col.names=F, row.names=F)
write.table(data.frame(x=row.names(as.data.frame(res_HHHHvsAAAA_f2)), y=-log10(as.data.frame(res_HHHHvsAAAA_f2)[,5])), 
    file="~/reciprocal_t/analysis/GO_enrich/dge_f2.txt",
    sep=",", quote=FALSE, col.names=F, row.names=F)
write.table(data.frame(x=row.names(as.data.frame(res_HHHHvsAAAA_f3)), y=-log10(as.data.frame(res_HHHHvsAAAA_f3)[,5])), 
    file="~/reciprocal_t/analysis/GO_enrich/dge_f3.txt",
    sep=",", quote=FALSE, col.names=F, row.names=F)

############################################
#
# SNP plot
#
############################################

all_out <- read.table("~/reciprocal_t/analysis/AF_change.txt", header=TRUE)

mean.out <- all_out[which(all_out$sig_all == TRUE),]

#looking at what more strict cutoff produces. pattern is consistent.
#f1 <- all_out[which(all_out$pval_f1 < quantile(all_out$pval_f1, 0.025, na.rm=TRUE)),]
#f2 <- all_out[which(all_out$pval_f2 < quantile(all_out$pval_f2, 0.025, na.rm=TRUE)),]
#f3 <- all_out[which(all_out$pval_f3 < quantile(all_out$pval_f3, 0.025, na.rm=TRUE)),]
#f1_f2 <- f1[which(f1$SNP %in% f2$SNP),]
#mean.out <- f1_f2[which(f1_f2$SNP %in% f3$SNP),]

pdf("~/reciprocal_t/figures/snp_AM.pdf", height = 3.75, width = 1.75)

par(mfrow = c(3, 1))
par(cex = 0.6)
par(mar = c(0, 3.7, 0, 0), oma = c(3.8, 0, 0.5, 0.9))
par(tcl = -0.25)
par(mgp = c(1.5, 0.3, 0))

# and for AA in HH
x <- mean.out$HHHH_F1_mean - mean.out$AAAA_F1_mean
y <- mean.out$HHAA_F1_mean - mean.out$AAAA_F1_mean
plot(x =x,
    y = y, 
    col=alpha("black", 0.4), bg=alpha("gray38", 1), pch=21, 
    ylim=c(-0.85,0.6), xlim=c(-.8, .8),
    xlab="",
    ylab="",
    cex.lab=1, cex.main=1.4, cex=1.2, xaxt='n', yaxt='n')
axis(2, at=c(-0.5, 0, 0.5), cex.axis=0.8)

nonAdapt_AA_f1 <- which((y > 0 & x < 0) | (y < 0 & x > 0))
points(x= x[nonAdapt_AA_f1], y=y[nonAdapt_AA_f1],
    col=alpha("black", 0.4), bg=alpha("gray86", 1), pch=21, cex=1.2)
Adapt_AA_f1 <- which((y < 0 & x < 0) | (y > 0 & x > 0))

abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="firebrick3", lwd=1.5)
#title(xlab="Adaptive divergence in allele frequency:\nGH - AM",mgp=c(3, 0.5, 0), cex.lab=1.3)
corval <- round(cor.test(x,y, method="pearson")$estimate, 2)
text(x=0.39, y=-0.6, bquote(rho == .(corval)), cex=1.2)

# F2
x <- mean.out$HHHH_F2_mean - mean.out$AAAA_F2_mean
y <- mean.out$HHAA_F2_mean - mean.out$AAAA_F2_mean
#summary(lm(y~x))

plot(x = x,
    y = y, 
    col=alpha("black", 0.4), bg=alpha("gray38", 1), pch=21, 
    ylim=c(-0.85,0.6), xlim=c(-.8, .8),
    xlab="",
    ylab="",
    cex.lab=1, cex.main=1.4, cex=1.2, xaxt='n', yaxt='n')
axis(2, at=c(-0.5, 0, 0.5), cex.axis=0.8)

nonAdapt_AA_f2 <- which((y > 0 & x < 0) | (y < 0 & x > 0))
points(x= x[nonAdapt_AA_f2], y=y[nonAdapt_AA_f2],
    col=alpha("black", 0.4), bg=alpha("gray86", 1), pch=21, cex=1.2)
Adapt_AA_f2 <- which((y < 0 & x < 0) | (y > 0 & x > 0))

abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="firebrick3", lwd=1.5)

corval <- round(cor.test(x,y, method="pearson")$estimate, 2)
text(x=0.39, y=-0.6, bquote(rho == .(corval)), cex=1.2)

 par(xpd=TRUE)
mtext(text = "Change in AM frequency due to treatment:\n(AM in GH) - AM",
      side = 2,#side 1 = bottom
      line = 1.7, cex=0.6)
 par(xpd=FALSE)

## F3
x <- mean.out$HHHH_F3_mean - mean.out$AAAA_F3_mean
y <- mean.out$HHAA_F3_mean - mean.out$AAAA_F3_mean
plot(x = x,
    y = y, 
    col=alpha("black", 0.4), bg=alpha("gray38", 1), pch=21, 
    ylim=c(-0.85,0.6), xlim=c(-.8, .8),
    xlab="",
    ylab="",
    cex.lab=1, cex.main=1.4, cex=1.2, xaxt='n', yaxt='n')
axis(2, at=c(-0.5, 0, 0.5), cex.axis=0.8)
axis(1, cex.axis=0.8)

nonAdapt_AA_f3 <- which((y > 0 & x < 0) | (y < 0 & x > 0))
points(x= x[nonAdapt_AA_f3], y=y[nonAdapt_AA_f3],
    col=alpha("black", 0.4), bg=alpha("gray86", 1), pch=21, cex=1.2)
Adapt_AA_f3 <- which((y < 0 & x < 0) | (y > 0 & x > 0))

abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="firebrick3", lwd=1.5)
#title(xlab="Adaptive divergence in allele frequency:\nGH - AM",mgp=c(3, 0.5, 0), cex.lab=1.3)
corval <- round(cor.test(x,y, method="pearson")$estimate, 2)
text(x=0.39, y=-0.6, bquote(rho == .(corval)), cex=1.2)



 par(xpd=TRUE)
mtext(text = "Adaptive divergence in allele\nfrequency: GH-AM",
      side = 1,#side 1 = bottom
      line = 2.3, cex=0.55)
 par(xpd=FALSE)

dev.off()



################
# HH
################

pdf("~/reciprocal_t/figures/snp_GH.pdf", height = 3.75, width = 1.75)

par(mfrow = c(3, 1))
par(cex = 0.6)
par(mar = c(0, 3.7, 0, 0), oma = c(3.8, 0, 0.5, 0.9))
par(tcl = -0.25)
par(mgp = c(1.5, 0.3, 0))

x <- mean.out$AAAA_F1_mean - mean.out$HHHH_F1_mean
y <- mean.out$AAHH_F1_mean - mean.out$HHHH_F1_mean
plot(x =x,
    y = y, 
    col=alpha("black", 0.4), bg=alpha("gray38", 1), pch=21, 
    ylim=c(-0.85,0.6), xlim=c(-.8, .8),
    xlab="",
    ylab="",
    cex.lab=1, cex.main=1.4, cex=1.2, xaxt='n', yaxt='n')
axis(2, at=c(-0.5, 0, 0.5), cex.axis=0.8)

nonAdapt_HH_f1 <- which((y > 0 & x < 0) | (y < 0 & x > 0))
points(x= x[nonAdapt_HH_f1], y=y[nonAdapt_HH_f1],
    col=alpha("black", 0.4), bg=alpha("gray86", 1), pch=21, cex=1.2)
Adapt_HH_f1 <- which((y < 0 & x < 0) | (y > 0 & x > 0))

abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="firebrick3", lwd=2)
#title(xlab="Adaptive divergence in allele frequency:\nAM - GH",mgp=c(3, 0.5, 0), cex.lab=1.3)
title(xlab="",mgp=c(3, 0.5, 0), cex.lab=1.3)
corval <- round(cor.test(x,y, method="pearson")$estimate, 2)
text(x=0.39, y=-0.6, bquote(rho == .(corval)), cex=1.2)


## F2
x <- mean.out$AAAA_F2_mean - mean.out$HHHH_F2_mean
y <- mean.out$AAHH_F2_mean - mean.out$HHHH_F2_mean
plot(x = x,
    y = y, 
    col=alpha("black", 0.4), bg=alpha("gray38", 1), pch=21, 
    ylim=c(-0.85,0.6), xlim=c(-.8, .8),
    xlab="",
    ylab="",
    cex.lab=1, cex.main=1.4, cex=1.2, xaxt='n', yaxt='n')
axis(2, at=c(-0.5, 0, 0.5), cex.axis=0.8)

nonAdapt_HH_f2 <- which((y > 0 & x < 0) | (y < 0 & x > 0))
points(x= x[nonAdapt_HH_f2], y=y[nonAdapt_HH_f2],
    col=alpha("black", 0.4), bg=alpha("gray86", 1), pch=21, cex=1.2)
Adapt_HH_f2 <- which((y < 0 & x < 0) | (y > 0 & x > 0))

abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="firebrick3", lwd=2)
#title(xlab="Adaptive divergence in allele frequency: AM - GH",mgp=c(3, 0.5, 0), cex.lab=1.5)
corval <- round(cor.test(x,y, method="pearson")$estimate, 2)
text(x=0.39, y=-0.6, bquote(rho == .(corval)), cex=1.2)

 par(xpd=TRUE)
mtext(text = "Change in GH frequency due to treatment:\n(GH in AM) - GH",
      side = 2,#side 1 = bottom
      line = 1.7, cex=0.6)
 par(xpd=FALSE)

## F3
x <- mean.out$AAAA_F3_mean - mean.out$HHHH_F3_mean
y <- mean.out$AAHH_F3_mean - mean.out$HHHH_F3_mean
plot(x = x,
    y = y, 
    col=alpha("black", 0.4), bg=alpha("gray38", 1), pch=21, 
    ylim=c(-0.85,0.6), xlim=c(-.8, .8),
    xlab="",
    ylab="",
    cex.lab=1, cex.main=1.4, cex=1.2, xaxt='n', yaxt='n')
axis(2, at=c(-0.5, 0, 0.5), cex.axis=0.8)
axis(1, cex.axis=0.8)

nonAdapt_HH_f3 <- which((y > 0 & x < 0) | (y < 0 & x > 0))
points(x= x[nonAdapt_HH_f3], y=y[nonAdapt_HH_f3],
    col=alpha("black", 0.4), bg=alpha("gray86", 1), pch=21, cex=1.2)
Adapt_HH_f3 <- which((y < 0 & x < 0) | (y > 0 & x > 0))

abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="firebrick3", lwd=1.5)
#title(xlab="Adaptive divergence in allele frequency:\nAM - GH",mgp=c(3, 0.5, 0), cex.lab=1.3)
corval <- round(cor.test(x,y, method="pearson")$estimate, 2)
text(x=0.39, y=-0.6, bquote(rho == .(corval)), cex=1.2)

 par(xpd=TRUE)
mtext(text = "Adaptive divergence in allele\nfrequency: AM-GH",
      side = 1,#side 1 = bottom
      line = 2.3, cex=0.55)
 par(xpd=FALSE)


dev.off()

# prop adaptive
length(Adapt_HH_f1)/sum(length(nonAdapt_HH_f1), length(Adapt_HH_f1))
# [1] 0.6859354
length(Adapt_HH_f2)/sum(length(nonAdapt_HH_f2), length(Adapt_HH_f2))
# [1] 0.7542563
length(Adapt_HH_f3)/sum(length(nonAdapt_HH_f3), length(Adapt_HH_f3))
# [1] 0.7613022
length(Adapt_AA_f1)/sum(length(nonAdapt_AA_f1), length(Adapt_AA_f1))
# [1] 0.6826923
length(Adapt_AA_f2)/sum(length(nonAdapt_AA_f2), length(Adapt_AA_f2))
# [1] 0.6409361
length(Adapt_AA_f3)/sum(length(nonAdapt_AA_f3), length(Adapt_AA_f3))
# [1] 0.6553224


prop.test(c(length(Adapt_HH_f1),length(Adapt_HH_f2),length(Adapt_HH_f3)),
          c(sum(length(nonAdapt_HH_f1), length(Adapt_HH_f1)),
            sum(length(nonAdapt_HH_f2),length(Adapt_HH_f2)),
            sum(length(nonAdapt_HH_f3), length(Adapt_HH_f3))),correct=FALSE)
# X-squared = 301.36, df = 2, p-value < 2.2e-16
prop.test(c(length(Adapt_HH_f1),length(Adapt_HH_f2)),
          c(sum(length(nonAdapt_HH_f1), length(Adapt_HH_f1)),
            sum(length(nonAdapt_HH_f2),length(Adapt_HH_f2))
            ),correct=FALSE)#
# X-squared = 196.63, df = 1, p-value < 2.2e-16
prop.test(c(length(Adapt_HH_f2),length(Adapt_HH_f3)),
          c(sum(length(nonAdapt_HH_f2), length(Adapt_HH_f2)),
            sum(length(nonAdapt_HH_f3),length(Adapt_HH_f3))
            ),correct=FALSE)
# X-squared = 2.2819, df = 1, p-value = 0.1309

prop.test(c(length(Adapt_HH_f1),length(Adapt_HH_f3)),
          c(sum(length(nonAdapt_HH_f1), length(Adapt_HH_f1)),
            sum(length(nonAdapt_HH_f3),length(Adapt_HH_f3))
            ),correct=FALSE)
# X-squared = 240.49, df = 1, p-value < 2.2e-16

# and for AA
prop.test(c(length(Adapt_AA_f1),length(Adapt_AA_f2),length(Adapt_AA_f3)),
          c(sum(length(nonAdapt_AA_f1), length(Adapt_AA_f1)),
            sum(length(nonAdapt_AA_f2), length(Adapt_AA_f2)),
            sum(length(nonAdapt_AA_f3), length(Adapt_AA_f3))),correct=FALSE)
# X-squared = 68.563, df = 2, p-value = 1.293e-15
prop.test(c(length(Adapt_AA_f1),length(Adapt_AA_f2)),
          c(sum(length(nonAdapt_AA_f1), length(Adapt_AA_f1)),
            sum(length(nonAdapt_AA_f2), length(Adapt_AA_f2))
            ),correct=FALSE)#
# X-squared = 66.665, df = 1, p-value = 3.218e-16
prop.test(c(length(Adapt_AA_f1),length(Adapt_AA_f3)),
          c(sum(length(nonAdapt_AA_f1), length(Adapt_AA_f1)),
            sum(length(nonAdapt_AA_f3), length(Adapt_AA_f3))
            ),correct=FALSE)#
# X-squared = 29.008, df = 1, p-value = 7.208e-08

prop.test(c(length(Adapt_AA_f2),length(Adapt_AA_f3)),
          c(sum(length(nonAdapt_AA_f2), length(Adapt_AA_f2)),
            sum(length(nonAdapt_AA_f3), length(Adapt_AA_f3))
            ),correct=FALSE)
# X-squared = 7.8106, df = 1, p-value = 0.005194

####################################
####################################
#
# compare cmh and dge changes
#
####################################
####################################


all_out$CHR <- sapply(strsplit(as.character(all_out$SNP), ":"), `[`, 1)

HHHHvsAAAA_f1 <- as.data.frame(res_HHHHvsAAAA_f1)
HHHHvsAAAA_f1$CHR <- row.names(HHHHvsAAAA_f1)
f1 <- merge(all_out, HHHHvsAAAA_f1, by="CHR", all=F)

f1$pvalue <- -log10(f1$pvalue)
f1$pval_f1 <- -log10(f1$pval_f1)
f1 <- (f1[which(!is.na(f1$pval_f1)),])

HHHHvsAAAA_f2 <- as.data.frame(res_HHHHvsAAAA_f2)
HHHHvsAAAA_f2$CHR <- row.names(HHHHvsAAAA_f2)
f2 <- merge(all_out, HHHHvsAAAA_f2, by="CHR")
f2$pvalue <- -log10(f2$pvalue)
f2$pval_f2 <- -log10(f2$pval_f2)
f2 <- (f2[which(!is.na(f2$pval_f2)),])
f2 <- (f2[which(is.finite(f2$pval_f2)),])

HHHHvsAAAA_f3 <- as.data.frame(res_HHHHvsAAAA_f3)
HHHHvsAAAA_f3$CHR <- row.names(HHHHvsAAAA_f3)
f3 <- merge(all_out, HHHHvsAAAA_f3, by="CHR")

f3$pvalue <- -log10(f3$pvalue)
f3$pval_f3 <- -log10(f3$pval_f3)
f3 <- (f3[which(!is.na(f3$pval_f3)),])
f3 <- (f3[which(is.finite(f3$pval_f3)),])

res_HHHHvsAAAA_f1
res_HHHHvsAAAA_f2
res_HHHHvsAAAA_f3

# look for the correlation between dge and af change. Does one drive the other?
png("~/reciprocal_t/figures/dge_af_cor.png", height = 2.5, width = 9, units="in", res=300)

par(mfrow = c(1, 3),mar = c(5, 5, 2, 3))

plot(x= (f1$pval_f1), y=(f1$pvalue),
    ylab = c("Gene expression: -log10(pvalue)"),
    xlab = c("Allele freq. change: -log10(pvalue)"),
    main = c("F1; R-sq = 0.001"),
    pch=19, col=alpha("black", 0.1))

summary(lm(f1$pvalue ~ f1$pval_f1))
# slope = 0.0042087
# rsq = 0.001024

plot(x= (f2$pval_f2), y=(f2$pvalue),
    ylab = c("Gene expression: -log10(pvalue)"),
    xlab = c("Allele freq. change: -log10(pvalue)"),
    main = c("F2; R-sq = 0.0001"),
    pch=19, col=alpha("black", 0.1))
summary(lm(f2$pvalue ~ f2$pval_f2, na.action = na.exclude))
# slope = 0.0010405
# rsq = 0.0001274

plot(x= f3$pval_f3, y=f3$pvalue,
    ylab = c("Gene expression: -log10(pvalue)"),
    xlab = c("Allele freq. change: -log10(pvalue)"),
    main = c("F3; R-sq = 0.002"),
    pch=19, col=alpha("black", 0.1))
summary(lm((f3$pvalue) ~ (f3$pval_f3)))
# slope = 0.0036867
# rsq = 0.001542
dev.off()


