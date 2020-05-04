

library(tximport)
library(tidyr)
library(DESeq2)
library(ggplot2)
library(ggpubr)


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
### 
###
###################################################################################
###################################################################################

files <- all_files

# make sure readr package is installed
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)

# making labels for sample table
id <- separate(data=samples, col=V1, sep="_", into = c("Population", "Generation", "Replicate"))
id$Treatment <- substr(id$Population, 1,2)
id$Line <- substr(id$Population, 3,4)
id$group <- paste(id$Treatment, id$Line, sep="")

sampleTable <- data.frame(
        sampleName = samples$V1,
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
keep <- apply(counts(dds), 1, function(x) {ifelse(length(which(x > 10)) > 43, TRUE, FALSE)})
dds <- dds[keep,]
nrow(dds)
#[1] 22545

# then transform
rld_f1<-vst(dds,blind=TRUE)


#######################################################
#######################################################
##
## PCA
##
#######################################################
#######################################################

#data <- plotPCA(rld_f1, intgroup=c("Treatment","Line"), returnData=TRUE, ntop = 22545)
#data <- plotPCA(rld_f1, intgroup=c("Treatment","Line"), returnData=TRUE, ntop = 500)

#rv <- rowVars(assay(rld_f1))
# select the ntop genes by variance
#select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
# perform a PCA on the data in assay(x) for the selected genes
#pca <- prcomp(t(assay(rld_f1)[select,]))
pca <- prcomp(t(assay(rld_f1)))

data <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3], PC4=pca$x[,4],
            Treatment = id$Treatment,
            Line = id$Line,
            name = colnames(rld_f1))

percentVar <- round(100* (pca$sdev^2)/sum(pca$sdev^2), digits=1)[1:10]

data$Treatment <- gsub("AA", "AM", data$Treatment)
data$Treatment <- gsub("HH", "GH", data$Treatment)
data$Line <- gsub("AA", "AM", data$Line)
data$Line <- gsub("HH", "GH", data$Line)
data$Generation <- substr(row.names(data),6,7)

f1in <- data[which(data$Generation == "F1"),]
f2in <- data[which(data$Generation == "F2"),]
f3in <- data[which(data$Generation == "F3"),]

a <- ggplot(f1in, aes(PC1, PC2, fill=Treatment, shape=Line, size=Line)) +
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
       #  stat_conf_ellipse(aes(color=Treatment, lty=Line), size=1, level=0.99,) +
       # scale_color_manual(values=c('cornflowerblue', 'brown2'))


b <- ggplot(f2in, aes(PC1, PC2, fill=Treatment, shape=Line, size=Line)) +
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

c <- ggplot(f3in, aes(PC1, PC2, fill=Treatment, shape=Line, size=Line)) +
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

png("~/reciprocal_t/figures/pca_RNA_combined_gens.png", res=300, height=4, width=11, units="in")

ggarrange(a, b,c, nrow = 1, ncol=3, common.legend=TRUE)

dev.off()

##### plot pc3 and pc 4

a2 <- ggplot(f1in, aes(PC3, PC4, fill=Treatment, shape=Line, size=Line)) +
        geom_point() +
        xlab(paste0("PC3: ",percentVar[3],"% variance")) +
        ylab(paste0("PC4: ",percentVar[4],"% variance")) +
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
       #  stat_conf_ellipse(aes(color=Treatment, lty=Line), size=1, level=0.99,) +
       # scale_color_manual(values=c('cornflowerblue', 'brown2'))


b2 <- ggplot(f2in, aes(PC3, PC4, fill=Treatment, shape=Line, size=Line)) +
        geom_point() +
        xlab(paste0("PC3: ",percentVar[3],"% variance")) +
        ylab(paste0("PC4: ",percentVar[4],"% variance")) +
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

c2 <- ggplot(f3in, aes(PC3, PC4, fill=Treatment, shape=Line, size=Line)) +
        geom_point() +
        xlab(paste0("PC3: ",percentVar[3],"% variance")) +
        ylab(paste0("PC4: ",percentVar[4],"% variance")) +
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

png("~/reciprocal_t/figures/pca_RNA_combined_gens_PC3_4.png", res=300, height=4, width=11, units="in")

ggarrange(a2, b2,c2, nrow = 1, ncol=3, common.legend=TRUE)

dev.off()

# write output for pca loadings for go enrichment:

loadings.out <- data.frame(
          SNP = rownames(pca$rotation), 
          PC1= pca$rotation[,1])
write.table(file="~/reciprocal_t/analysis/GO_enrich/dge_contrib_PC1.txt", loadings.out, sep="\t", 
              row.names=FALSE, quote=FALSE)

loadings.out <- data.frame(
          SNP = rownames(pca$rotation), 
          PC2= pca$rotation[,2])
write.table(file="~/reciprocal_t/analysis/GO_enrich/dge_contrib_PC2.txt", loadings.out, sep="\t", 
              row.names=FALSE, quote=FALSE)

#############################################
#############################################
##### snps
#############################################
#############################################

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

# all gens together
pcaResult <- prcomp(freqs)

percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

data <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3], PC4=pca$x[,4],
            Treatment = id$Treatment,
            Line = id$Line,
            name = colnames(rld_f1))

percentVar <- round(100* (pcaResult$sdev^2)/sum(pcaResult$sdev^2), digits=1)[1:10]

data <- data.frame(id=pops, Line=substr(pops, 3,4), 
    Treatment=substr(pops, 1,2),
        Generations=substr(pops, 6,7),
        PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,2],
        PC3 = pcaResult$x[,3],  PC4= pcaResult$x[,4])


dat_f1 <- data[which(data$Generation == "F1"),]
dat_f2 <- data[which(data$Generation == "F2"),]
dat_f3 <- data[which(data$Generation == "F3"),]


d <- ggplot(dat_f1, aes(PC1, PC2, fill=Treatment, shape=Line, size=Line)) +
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

e <- ggplot(dat_f2, aes(PC1, PC2, fill=Treatment, shape=Line, size=Line)) +
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

f <- ggplot(dat_f3, aes(PC1, PC2, fill=Treatment, shape=Line, size=Line)) +
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

ggarrange(d, e, f, nrow = 1, ncol=3, common.legend=TRUE)

dev.off()

pdf("~/reciprocal_t/figures/fig1_pca.pdf",height=8, width=11)

ggarrange(d,e,f,a,b,c, nrow = 2, ncol=3, common.legend=TRUE, labels="AUTO", legend="bottom")

dev.off()


### pc 3 and 4


d2 <- ggplot(dat_f1, aes(PC3, PC4, fill=Treatment, shape=Line, size=Line)) +
        geom_point() +
        xlab(paste0("PC3: ",percentVar[3],"% variance")) +
        ylab(paste0("PC4: ",percentVar[4],"% variance")) +
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

e2 <- ggplot(dat_f2, aes(PC3, PC4, fill=Treatment, shape=Line, size=Line)) +
        geom_point() +
        xlab(paste0("PC3: ",percentVar[3],"% variance")) +
        ylab(paste0("PC4: ",percentVar[4],"% variance")) +
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

f2 <- ggplot(dat_f3, aes(PC3, PC4, fill=Treatment, shape=Line, size=Line)) +
        geom_point() +
        xlab(paste0("PC3: ",percentVar[3],"% variance")) +
        ylab(paste0("PC4: ",percentVar[4],"% variance")) +
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


png("~/reciprocal_t/figures/pca_SNP_PC3_4.png", res=300, height=4, width=11, units="in")

ggarrange(d2, e2, f2, nrow = 1, ncol=3, common.legend=TRUE)

dev.off()

