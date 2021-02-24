


# run deseq2 to get plastic, genetic, etc. changes in gene expression.

#############################################
#############################################
#### F1
#############################################
#############################################

library(tximport)
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(tximport)
# running gene expression data first here:

#locate the directory containing the files.
dir <- "~/reciprocal_t/analysis/salmon"
samples <- as.vector(list.files(dir))

# read in table with sample ids
# samples <- read.table("~/tonsa_genomics/analysis/sample_id.txt", header=FALSE)

samples_in <- samples[grep("F1", samples, invert=F)]
#samples_in <- samples_in[grep("AAHH", samples_in, invert=T)]
#samples_in <- samples_in[grep("HHAA", samples_in, invert=T)]
samples <- as.data.frame(samples_in)
colnames(samples) <- c("V1")
# now point to quant files
all_files <- file.path(dir, samples$V1, "quant.sf")
names(all_files) <- samples$V1

# associate transcripts with gene ids
# make a data.frame called tx2gene with two columns: 1) transcript ID and 2) gene ID. col order is important
gene_tran <- read.table("/data/copepods/tonsa_transcriptome/Atonsa_gen_trans_agp_gff/Atonsa_transcript_to_gene", header=FALSE)
# then convert
tx2gene <- data.frame(transcript=gene_tran$V2, gene=gene_tran$V1)


# make sure readr package is installed
txi <- tximport(all_files, type = "salmon", tx2gene = tx2gene)
names(txi)

# making labels for sample table

id <- separate(data=samples, col=V1, sep="_", into = c("Population","generation", "Replicate"))
id$Line <- substr(id$Population, 3,4)
id$group <- paste(id$Treatment, id$Line, sep="")

sampleTable <- data.frame(
        sampleName = samples$V1,
        Line = id$Line,
        Population = id$Population,
        Replicate = id$Replicate)

# assign row names to sample table
rownames(sampleTable) <- colnames(txi$counts)

# convert to DEseq
dds <- DESeqDataSetFromTximport(txi,
                colData=sampleTable,
                design = ~ Population)
                      #Treatment:Line + Treatment*Line*Generation)

# remove rows where count < 10 in more than 75% of samples
keep <- apply(counts(dds), 1, function(x) {ifelse(length(which(x > 10)) > (nrow(samples)*0.75), TRUE, FALSE)})
dds <- dds[keep,]
nrow(dds)
#[1] 24362

dds <- estimateSizeFactors(dds)
# save counts:
write.table(counts(dds, normalized=TRUE), file="~/reciprocal_t/analysis/DE_counts_F1.txt", sep="\t", quote=FALSE)

# then rolog transform
rld<-rlog(dds,blind=TRUE)

dat_rld=as.data.frame(assay(rld))
colnames(dat_rld)<-colnames(dds)

#write.table(dat_rld, file="~/tonsa_genomics/analysis/rlog_counts.txt", sep="\t", quote=FALSE)

##########################
#### run DEseq
##########################

# change levels. so I know what the reference level is.
dds$Population <- relevel(dds$Population, ref = "AAAA")

dds <- DESeq(dds, parallel=T)
res <- results(dds, alpha=0.1)

dds <- estimateSizeFactors(dds)

write.table(counts(dds, normalized=TRUE), file="~/reciprocal_t/analysis/norm_counts_F1.txt", sep="\t", quote=FALSE)


res_summary <- results(dds, contrast=c("Population", "HHHH","AAAA"))
summary(res_summary, alpha=0.05)

write.table(as.data.frame(res_summary), file="~/reciprocal_t/analysis/DGE_AM_vs_GH_results_F1.txt", sep="\t", quote=FALSE)



#############################################
## plasticity only
#############################################

res_summary <- results(dds, contrast=c("Population", "HHAA","AAAA"))
summary(res_summary, alpha=0.05)

write.table(as.data.frame(res_summary), file="~/reciprocal_t/analysis/DGE_AM_plasticity_results_F1.txt", sep="\t", quote=FALSE)

###### GH plasticity

res_summary <- results(dds, contrast=c("Population", "HHHH", "AAHH"))
summary(res_summary, alpha=0.05)

write.table(as.data.frame(res_summary), file="~/reciprocal_t/analysis/DGE_GH_plasticity_results_F1.txt", sep="\t", quote=FALSE)

#############################################
## evolved
#############################################

res_summary <- results(dds, contrast=c("Population", "AAHH","AAAA"))
summary(res_summary, alpha=0.05)

write.table(as.data.frame(res_summary), file="~/reciprocal_t/analysis/DGE_AM_GH_in_AM_results_F1.txt", sep="\t", quote=FALSE)


###### GH plasticity

res_summary <- results(dds, contrast=c("Population", "HHAA", "HHHH"))
summary(res_summary, alpha=0.05)

write.table(as.data.frame(res_summary), file="~/reciprocal_t/analysis/DGE_AM_GH_in_GH_results_F1.txt", sep="\t", quote=FALSE)

#############################################
#############################################
################ F2
#############################################
#############################################

#locate the directory containing the files.
dir <- "~/reciprocal_t/analysis/salmon"
samples <- as.vector(list.files(dir))

# read in table with sample ids
# samples <- read.table("~/tonsa_genomics/analysis/sample_id.txt", header=FALSE)

samples_in <- samples[grep("F2", samples, invert=F)]

samples <- as.data.frame(samples_in)
colnames(samples) <- c("V1")
# now point to quant files
all_files <- file.path(dir, samples$V1, "quant.sf")
names(all_files) <- samples$V1

# associate transcripts with gene ids
# make a data.frame called tx2gene with two columns: 1) transcript ID and 2) gene ID. col order is important
gene_tran <- read.table("/data/copepods/tonsa_transcriptome/Atonsa_gen_trans_agp_gff/Atonsa_transcript_to_gene", header=FALSE)
# then convert
tx2gene <- data.frame(transcript=gene_tran$V2, gene=gene_tran$V1)


# make sure readr package is installed
txi <- tximport(all_files, type = "salmon", tx2gene = tx2gene)
names(txi)

# making labels for sample table

id <- separate(data=samples, col=V1, sep="_", into = c("Population","generation", "Replicate"))
id$Line <- substr(id$Population, 3,4)
id$group <- paste(id$Treatment, id$Line, sep="")

sampleTable <- data.frame(
        sampleName = samples$V1,
        Line = id$Line,
        Population = id$Population,
        Replicate = id$Replicate)

# assign row names to sample table
rownames(sampleTable) <- colnames(txi$counts)

# convert to DEseq
dds <- DESeqDataSetFromTximport(txi,
                colData=sampleTable,
                design = ~ Population)
                      #Treatment:Line + Treatment*Line*Generation)

# remove rows where count < 10 in more than 75% of samples
keep <- apply(counts(dds), 1, function(x) {ifelse(length(which(x > 10)) > (nrow(samples)*0.75), TRUE, FALSE)})
dds <- dds[keep,]
nrow(dds)
#[1] 25948
dds <- estimateSizeFactors(dds)
# save counts:
write.table(counts(dds, normalized=TRUE), file="~/reciprocal_t/analysis/DE_counts_F2.txt", sep="\t", quote=FALSE)

# then rolog transform
rld<-rlog(dds,blind=TRUE)

dat_rld=as.data.frame(assay(rld))
colnames(dat_rld)<-colnames(dds)


##########################
#### run DEseq
##########################

# change levels. so I know what the reference level is.
dds$Line <- relevel(dds$Population, ref = "AAAA")

dds <- DESeq(dds, parallel=T)
res <- results(dds, alpha=0.1)

dds <- estimateSizeFactors(dds)

write.table(counts(dds, normalized=TRUE), file="~/reciprocal_t/analysis/norm_counts_F2.txt", sep="\t", quote=FALSE)

res_summary <- results(dds, contrast=c("Population", "HHHH","AAAA"))
summary(res_summary, alpha=0.05)
write.table(as.data.frame(res_summary), file="~/reciprocal_t/analysis/DGE_AM_vs_GH_results_F2.txt", sep="\t", quote=FALSE)

#############################################
## plasticity
#############################################

res_summary <- results(dds, contrast=c("Population", "HHAA","AAAA"))
summary(res_summary, alpha=0.05)

write.table(as.data.frame(res_summary), file="~/reciprocal_t/analysis/DGE_AM_plasticity_results_F2.txt", sep="\t", quote=FALSE)

###### GH plasticity
res_summary <- results(dds, contrast=c("Population", "HHHH", "AAHH"))
summary(res_summary, alpha=0.05)

write.table(as.data.frame(res_summary), file="~/reciprocal_t/analysis/DGE_GH_plasticity_results_F2.txt", sep="\t", quote=FALSE)

#############################################
## evolved
#############################################

res_summary <- results(dds, contrast=c("Population", "AAHH","AAAA"))
summary(res_summary, alpha=0.05)

write.table(as.data.frame(res_summary), file="~/reciprocal_t/analysis/DGE_AM_GH_in_AM_results_F2.txt", sep="\t", quote=FALSE)

###### GH plasticity
res_summary <- results(dds, contrast=c("Population", "HHAA", "HHHH"))
summary(res_summary, alpha=0.05)

write.table(as.data.frame(res_summary), file="~/reciprocal_t/analysis/DGE_AM_GH_in_GH_results_F2.txt", sep="\t", quote=FALSE)




#############################################
#############################################
## F3
#############################################
#############################################

#locate the directory containing the files.
dir <- "~/reciprocal_t/analysis/salmon"
samples <- as.vector(list.files(dir))

# read in table with sample ids
# samples <- read.table("~/tonsa_genomics/analysis/sample_id.txt", header=FALSE)

samples_in <- samples[grep("F3", samples, invert=F)]

samples <- as.data.frame(samples_in)
colnames(samples) <- c("V1")
# now point to quant files
all_files <- file.path(dir, samples$V1, "quant.sf")
names(all_files) <- samples$V1

# associate transcripts with gene ids
# make a data.frame called tx2gene with two columns: 1) transcript ID and 2) gene ID. col order is important
gene_tran <- read.table("/data/copepods/tonsa_transcriptome/Atonsa_gen_trans_agp_gff/Atonsa_transcript_to_gene", header=FALSE)
# then convert
tx2gene <- data.frame(transcript=gene_tran$V2, gene=gene_tran$V1)


# make sure readr package is installed
txi <- tximport(all_files, type = "salmon", tx2gene = tx2gene)
names(txi)

# making labels for sample table

id <- separate(data=samples, col=V1, sep="_", into = c("Population","generation", "Replicate"))
id$Line <- substr(id$Population, 3,4)
id$group <- paste(id$Treatment, id$Line, sep="")

sampleTable <- data.frame(
        sampleName = samples$V1,
        Line = id$Line,
        Population = id$Population,
        Replicate = id$Replicate)

# assign row names to sample table
rownames(sampleTable) <- colnames(txi$counts)

# convert to DEseq
dds <- DESeqDataSetFromTximport(txi,
                colData=sampleTable,
                design = ~ Population)
                      #Treatment:Line + Treatment*Line*Generation)

# remove rows where count < 10 in more than 75% of samples
keep <- apply(counts(dds), 1, function(x) {ifelse(length(which(x > 10)) > (nrow(samples)*0.75), TRUE, FALSE)})
dds <- dds[keep,]
nrow(dds)
#[1] 25279
dds <- estimateSizeFactors(dds)
# save counts:
write.table(counts(dds, normalized=TRUE), file="~/reciprocal_t/analysis/DE_counts_F3.txt", sep="\t", quote=FALSE)

# then rolog transform
rld<-rlog(dds,blind=TRUE)

dat_rld=as.data.frame(assay(rld))
colnames(dat_rld)<-colnames(dds)

#write.table(dat_rld, file="~/tonsa_genomics/analysis/rlog_counts.txt", sep="\t", quote=FALSE)

##########################
#### run DEseq
##########################

# change levels. so I know what the reference level is.
dds$Line <- relevel(dds$Line, ref = "AAAA")

dds <- DESeq(dds, parallel=T)
res <- results(dds, alpha=0.1)

dds <- estimateSizeFactors(dds)

write.table(counts(dds, normalized=TRUE), file="~/reciprocal_t/analysis/norm_counts_F3.txt", sep="\t", quote=FALSE)


res_summary <- results(dds, contrast=c("Population", "HHHH","AAAA"))
summary(res_summary, alpha=0.05)

write.table(as.data.frame(res_summary), file="~/reciprocal_t/analysis/DGE_AM_vs_GH_results_F3.txt", sep="\t", quote=FALSE)

#############################################
## plasticity AM only
#############################################
res_summary <- results(dds, contrast=c("Population", "HHAA","AAAA"))
summary(res_summary, alpha=0.05)

write.table(as.data.frame(res_summary), file="~/reciprocal_t/analysis/DGE_AM_plasticity_results_F3.txt", sep="\t", quote=FALSE)

###### GH plasticity
res_summary <- results(dds, contrast=c("Population", "HHHH", "AAHH"))
summary(res_summary, alpha=0.05)

write.table(as.data.frame(res_summary), file="~/reciprocal_t/analysis/DGE_GH_plasticity_results_F3.txt", sep="\t", quote=FALSE)


#############################################
## evolved
#############################################
res_summary <- results(dds, contrast=c("Population", "AAHH","AAAA"))
summary(res_summary, alpha=0.05)

write.table(as.data.frame(res_summary), file="~/reciprocal_t/analysis/DGE_AM_GH_in_AM_results_F3.txt", sep="\t", quote=FALSE)

###### GH plasticity
res_summary <- results(dds, contrast=c("Population", "HHAA", "HHHH"))
summary(res_summary, alpha=0.05)

write.table(as.data.frame(res_summary), file="~/reciprocal_t/analysis/DGE_AM_GH_in_GH_results_F3.txt", sep="\t", quote=FALSE)








####################################################################################################
####################################################################################################
####################################################################################################
# plasticity plots
####################################################################################################
####################################################################################################
####################################################################################################

library(ggplot2)

# read in dge
dge_am <- read.csv("~/reciprocal_t/analysis/DGE_AM_plasticity_results_F1.txt", header=T, sep="\t")
colnames(dge_am) <- c("baseMean", colnames(dge_am)[2:ncol(dge_am)] )
dge_am$gene <- row.names(dge_am)

# gh plasticity
dge_gh <- read.csv("~/reciprocal_t/analysis/DGE_GH_plasticity_results_F1.txt", header=T, sep="\t")
colnames(dge_gh) <- c("baseMean", colnames(dge_gh)[2:ncol(dge_gh)] )
dge_gh$gene <- row.names(dge_gh)

names(dge_am)[names(dge_am) == 'log2FoldChange'] <- 'log2FoldChange.AM'
names(dge_gh)[names(dge_gh) == 'log2FoldChange'] <- 'log2FoldChange.GH'

names(dge_am)[names(dge_am) == 'padj'] <- 'padj.AM'
names(dge_gh)[names(dge_gh) == 'padj'] <- 'padj.GH'

# am vs gh at home:
dge_evol <- read.csv("~/reciprocal_t/analysis/DGE_AM_vs_GH_results_F1.txt", header=T, sep="\t")
colnames(dge_evol) <- c("baseMean", colnames(dge_evol)[2:ncol(dge_evol)] )
dge_evol$gene <- row.names(dge_evol)
names(dge_evol)[names(dge_evol) == 'log2FoldChange'] <- 'log2FoldChange.AM_vs_GH'

dge1 <- merge(dge_am, dge_gh, by="gene")

##########
# F2
##########
dge_am <- read.csv("~/reciprocal_t/analysis/DGE_AM_plasticity_results_F2.txt", header=T, sep="\t")
colnames(dge_am) <- c("baseMean", colnames(dge_am)[2:ncol(dge_am)] )
dge_am$gene <- row.names(dge_am)

# gh plasticity
dge_gh <- read.csv("~/reciprocal_t/analysis/DGE_GH_plasticity_results_F2.txt", header=T, sep="\t")
colnames(dge_gh) <- c("baseMean", colnames(dge_gh)[2:ncol(dge_gh)] )
dge_gh$gene <- row.names(dge_gh)

names(dge_am)[names(dge_am) == 'log2FoldChange'] <- 'log2FoldChange.AM'
names(dge_gh)[names(dge_gh) == 'log2FoldChange'] <- 'log2FoldChange.GH'

names(dge_am)[names(dge_am) == 'padj'] <- 'padj.AM'
names(dge_gh)[names(dge_gh) == 'padj'] <- 'padj.GH'

######################
# global loss in plasticity

dge2 <- merge(dge_gh, dge_am, by="gene")

# F3

dge_am <- read.csv("~/reciprocal_t/analysis/DGE_AM_plasticity_results_F3.txt", header=T, sep="\t")
colnames(dge_am) <- c("baseMean", colnames(dge_am)[2:ncol(dge_am)] )
dge_am$gene <- row.names(dge_am)

# gh plasticity
dge_gh <- read.csv("~/reciprocal_t/analysis/DGE_GH_plasticity_results_F3.txt", header=T, sep="\t")
colnames(dge_gh) <- c("baseMean", colnames(dge_gh)[2:ncol(dge_gh)] )
dge_gh$gene <- row.names(dge_gh)

names(dge_am)[names(dge_am) == 'log2FoldChange'] <- 'log2FoldChange.AM'
names(dge_gh)[names(dge_gh) == 'log2FoldChange'] <- 'log2FoldChange.GH'

names(dge_am)[names(dge_am) == 'padj'] <- 'padj.AM'
names(dge_gh)[names(dge_gh) == 'padj'] <- 'padj.GH'
######################
# global loss in plasticity

dge3 <- merge(dge_gh, dge_am, by="gene")

dge1$generation <- "F1"
dge2$generation <- "F2"
dge3$generation <- "F3"

dge <- rbind(rbind(dge1, dge2),dge3)

p <- ggplot(dge, aes(x=log2FoldChange.AM, y=log2FoldChange.GH)) +
  geom_point(size=2, shape=23) +
  xlab("AM plasticity: log2(AM in GH / AM in AM)") +
  ylab("GH plasticity: log2(GH in GH / Gh in AM)") +
   geom_smooth(method="lm") +
   geom_abline(intercept=0, slope=1) +
   facet_wrap(~generation)

ggsave(file = "~/reciprocal_t/figures/plasticity_generations.png", p, w=11.8, h=5.74)
#11.8 x 5.74


#####
# significant only
#####

dge_sig <- dge[which(dge$padj.AM < 0.05 | dge$padj.GH < 0.05),]

dge_sig$sig_class <- NA
dge_sig$sig_class[which(dge_sig$padj.AM < 0.05 & dge_sig$padj.GH >= 0.05)] <- "AM"
dge_sig$sig_class[which(dge_sig$padj.AM >= 0.05 & dge_sig$padj.GH < 0.05)] <- "GH"
dge_sig$sig_class[which(dge_sig$padj.AM < 0.05 & dge_sig$padj.GH < 0.05)] <- "both"

dge_sig <- dge_sig[which(!is.na(dge_sig$sig_class)),]
nrow(dge_sig)

p <- ggplot(dge_sig, aes(x=log2FoldChange.AM, y=log2FoldChange.GH, fill=sig_class,color=sig_class, shape=sig_class)) +
  geom_hline(yintercept=0, col="grey50") +
  geom_vline(xintercept=0, col="grey50") +
  geom_point(size=2, alpha=1, color="black") +
  scale_shape_manual(values = c(21, 22, 24)) +
  scale_color_manual(values=c("#6699CC","chartreuse4","#CC3333"))+
  scale_fill_manual(values=c("#6699CC","chartreuse4","#CC3333"))+
  xlab("AM plasticity: log2(AM in GH / AM in AM)") +
  ylab("GH plasticity: log2(GH in GH / Gh in AM)") +
   geom_abline(intercept=0, slope=1, show.legend = FALSE, lty=2) +
   #geom_smooth(method="lm", show.legend = FALSE) +
   geom_smooth(method="lm",aes(group=1), color="black", show.legend = FALSE) +
   facet_grid(rows=vars(generation)) +
   theme_bw() +
   guides(shape = guide_legend(override.aes = list(size = 5, alpha=1),
          fill = guide_legend(override.aes = list(size = 5, alpha=1))))+
theme(
  strip.background = element_blank(),
  strip.text = element_blank(),
  legend.position = "none"
)
#p

ggsave(file = "~/reciprocal_t/figures/plasticity_generations_sig_only.pdf", p, w=2.5, h=6.3)

# get slope, etc:

s_f1 <- dge_sig[which(dge_sig$generation == "F1"),]
s_f2 <- dge_sig[which(dge_sig$generation == "F2"),]
s_f3 <- dge_sig[which(dge_sig$generation == "F3"),]

summary(lm(s_f1$log2FoldChange.GH ~ s_f1$log2FoldChange.AM))
summary(lm(s_f2$log2FoldChange.GH ~ s_f2$log2FoldChange.AM))
summary(lm(s_f3$log2FoldChange.GH ~ s_f3$log2FoldChange.AM))








####################################################
####################################################
####################################################
# Ho et al analysis, sci advances 2020
####################################################
####################################################
####################################################

# need the following categories.
# 1. plastic within AM
# 2. plastic within GH
# 3. sig different between AMgh and GHgh
# 4. sig different between GHam and AMam

# 1. DGE_AM_plasticity_results_F1
# 2. DGE_GH_plasticity_results_F1
# 3. DGE_AM_GH_in_GH_results_F1
# 4. DGE_AM_GH_in_AM_results_F1

# PO forward = 1 but not 3
# PO reverse = 2 but not 4
# GN forward = 3
# GN reverse = 4

f1 <- read.csv("~/reciprocal_t/analysis/DGE_AM_vs_GH_results_F1.txt", header=T, sep="\t")
f2 <- read.csv("~/reciprocal_t/analysis/DGE_AM_vs_GH_results_F2.txt", header=T, sep="\t")
f3 <- read.csv("~/reciprocal_t/analysis/DGE_AM_vs_GH_results_F3.txt", header=T, sep="\t")

f1_sig <- row.names(f1)[which(f1$padj < 0.05)]
f2_sig <- row.names(f2)[which(f2$padj < 0.05)]
f3_sig <- row.names(f3)[which(f3$padj < 0.05)]

length(f1_sig)
length(f2_sig)
length(f3_sig)

### F1

# cat 1:
dge_am <- read.csv("~/reciprocal_t/analysis/DGE_AM_plasticity_results_F1.txt", header=T, sep="\t")
colnames(dge_am) <- c("baseMean", colnames(dge_am)[2:ncol(dge_am)] )
dge_am$gene <- row.names(dge_am)

dge_am <- dge_am[dge_am$gene %in% f1_sig,]

# cat 2: gh plasticity
dge_gh <- read.csv("~/reciprocal_t/analysis/DGE_GH_plasticity_results_F1.txt", header=T, sep="\t")
colnames(dge_gh) <- c("baseMean", colnames(dge_gh)[2:ncol(dge_gh)] )
dge_gh$gene <- row.names(dge_gh)

dge_gh <- dge_gh[dge_gh$gene %in% f1_sig,]

# cat 3
dge_gh_am_in_GH <- read.csv("~/reciprocal_t/analysis/DGE_AM_GH_in_GH_results_F1.txt", header=T, sep="\t")
colnames(dge_gh_am_in_GH) <- c("baseMean", colnames(dge_gh_am_in_GH)[2:ncol(dge_gh_am_in_GH)] )
dge_gh_am_in_GH$gene <- row.names(dge_gh_am_in_GH)
dge_gh_am_in_GH <- dge_gh_am_in_GH[dge_gh_am_in_GH$gene %in% f1_sig,]

# cat 4
dge_gh_am_in_AM <- read.csv("~/reciprocal_t/analysis/DGE_AM_GH_in_AM_results_F1.txt", header=T, sep="\t")
colnames(dge_gh_am_in_AM) <- c("baseMean", colnames(dge_gh_am_in_AM)[2:ncol(dge_gh_am_in_AM)] )
dge_gh_am_in_AM$gene <- row.names(dge_gh_am_in_AM)
dge_gh_am_in_AM <- dge_gh_am_in_AM[dge_gh_am_in_AM$gene %in% f1_sig,]

#######
## identify significant genes in each
dge_am_sig <- row.names(dge_am)[which(dge_am$padj < 0.05)]
dge_gh_sig <- row.names(dge_gh)[which(dge_gh$padj < 0.05)]


# evolved diff
forward_GN_f1 <- row.names(dge_gh_am_in_GH)[which(dge_gh_am_in_GH$padj < 0.05)]
reverse_GN_f1 <- row.names(dge_gh_am_in_AM)[which(dge_gh_am_in_AM$padj < 0.05)]

# remove the evolved from the plastic
forward_PO_f1 <- dge_am_sig[!dge_am_sig %in% forward_GN_f1]
reverse_PO_f1 <- dge_gh_sig[!dge_gh_sig %in% reverse_GN_f1]

length(forward_GN_f1)
length(reverse_GN_f1)
length(forward_PO_f1)
length(reverse_PO_f1)

forward_f1 <- length(forward_PO_f1)/length(forward_GN_f1)
reverse_f1 <- length(reverse_PO_f1)/length(reverse_GN_f1)

length(unique(c(forward_PO_f1, forward_GN_f1, reverse_GN_f1, reverse_PO_f1)))
# 87 % of the genes were accounted for.

######################
### F2
######################
# cat 1:
dge_am <- read.csv("~/reciprocal_t/analysis/DGE_AM_plasticity_results_F2.txt", header=T, sep="\t")
colnames(dge_am) <- c("baseMean", colnames(dge_am)[2:ncol(dge_am)] )
dge_am$gene <- row.names(dge_am)

dge_am <- dge_am[dge_am$gene %in% f2_sig,]

# cat 2: gh plasticity
dge_gh <- read.csv("~/reciprocal_t/analysis/DGE_GH_plasticity_results_F2.txt", header=T, sep="\t")
colnames(dge_gh) <- c("baseMean", colnames(dge_gh)[2:ncol(dge_gh)] )
dge_gh$gene <- row.names(dge_gh)

dge_gh <- dge_gh[dge_gh$gene %in% f2_sig,]

# cat 3
dge_gh_am_in_GH <- read.csv("~/reciprocal_t/analysis/DGE_AM_GH_in_GH_results_F2.txt", header=T, sep="\t")
colnames(dge_gh_am_in_GH) <- c("baseMean", colnames(dge_gh_am_in_GH)[2:ncol(dge_gh_am_in_GH)] )
dge_gh_am_in_GH$gene <- row.names(dge_gh_am_in_GH)
dge_gh_am_in_GH <- dge_gh_am_in_GH[dge_gh_am_in_GH$gene %in% f2_sig,]

# cat 4
dge_gh_am_in_AM <- read.csv("~/reciprocal_t/analysis/DGE_AM_GH_in_AM_results_F2.txt", header=T, sep="\t")
colnames(dge_gh_am_in_AM) <- c("baseMean", colnames(dge_gh_am_in_AM)[2:ncol(dge_gh_am_in_AM)] )
dge_gh_am_in_AM$gene <- row.names(dge_gh_am_in_AM)
dge_gh_am_in_AM <- dge_gh_am_in_AM[dge_gh_am_in_AM$gene %in% f2_sig,]

#######
## identify significant genes ine ach
dge_am_sig <- row.names(dge_am)[which(dge_am$padj < 0.05)]
dge_gh_sig <- row.names(dge_gh)[which(dge_gh$padj < 0.05)]

# evolved diff
forward_GN_f2 <- row.names(dge_gh_am_in_GH)[which(dge_gh_am_in_GH$padj < 0.05)]
reverse_GN_f2 <- row.names(dge_gh_am_in_AM)[which(dge_gh_am_in_AM$padj < 0.05)]

# remove the evolved from the plastic
forward_PO_f2 <- dge_am_sig[!dge_am_sig %in% forward_GN_f2]
reverse_PO_f2 <- dge_gh_sig[!dge_gh_sig %in% reverse_GN_f2]

length(forward_GN_f2)
length(reverse_GN_f2)
length(forward_PO_f2)
length(reverse_PO_f2)

forward_f2 <- length(forward_PO_f2)/length(forward_GN_f2)
reverse_f2 <- length(reverse_PO_f2)/length(reverse_GN_f2)

######################
### F3
######################

# cat 1:
dge_am <- read.csv("~/reciprocal_t/analysis/DGE_AM_plasticity_results_F3.txt", header=T, sep="\t")
colnames(dge_am) <- c("baseMean", colnames(dge_am)[2:ncol(dge_am)] )
dge_am$gene <- row.names(dge_am)

dge_am <- dge_am[dge_am$gene %in% f3_sig,]

# cat 2: gh plasticity
dge_gh <- read.csv("~/reciprocal_t/analysis/DGE_GH_plasticity_results_F3.txt", header=T, sep="\t")
colnames(dge_gh) <- c("baseMean", colnames(dge_gh)[2:ncol(dge_gh)] )
dge_gh$gene <- row.names(dge_gh)

dge_gh <- dge_gh[dge_gh$gene %in% f3_sig,]

# cat 3
dge_gh_am_in_GH <- read.csv("~/reciprocal_t/analysis/DGE_AM_GH_in_GH_results_F3.txt", header=T, sep="\t")
colnames(dge_gh_am_in_GH) <- c("baseMean", colnames(dge_gh_am_in_GH)[2:ncol(dge_gh_am_in_GH)] )
dge_gh_am_in_GH$gene <- row.names(dge_gh_am_in_GH)
dge_gh_am_in_GH <- dge_gh_am_in_GH[dge_gh_am_in_GH$gene %in% f3_sig,]

# cat 4
dge_gh_am_in_AM <- read.csv("~/reciprocal_t/analysis/DGE_AM_GH_in_AM_results_F3.txt", header=T, sep="\t")
colnames(dge_gh_am_in_AM) <- c("baseMean", colnames(dge_gh_am_in_AM)[2:ncol(dge_gh_am_in_AM)] )
dge_gh_am_in_AM$gene <- row.names(dge_gh_am_in_AM)
dge_gh_am_in_AM <- dge_gh_am_in_AM[dge_gh_am_in_AM$gene %in% f3_sig,]

#######
## identify significant genes ine ach
dge_am_sig <- row.names(dge_am)[which(dge_am$padj < 0.05)]
dge_gh_sig <- row.names(dge_gh)[which(dge_gh$padj < 0.05)]

# evolved diff
forward_GN_f3 <- row.names(dge_gh_am_in_GH)[which(dge_gh_am_in_GH$padj < 0.05)]
reverse_GN_f3 <- row.names(dge_gh_am_in_AM)[which(dge_gh_am_in_AM$padj < 0.05)]

# remove the evolved from the plastic
forward_PO_f3 <- dge_am_sig[!dge_am_sig %in% forward_GN_f3]
reverse_PO_f3 <- dge_gh_sig[!dge_gh_sig %in% reverse_GN_f3]

length(forward_GN_f3)
length(reverse_GN_f3)
length(forward_PO_f3)
length(reverse_PO_f3)

forward_f3 <- length(forward_PO_f3)/length(forward_GN_f3)
reverse_f3 <- length(reverse_PO_f3)/length(reverse_GN_f3)

pdat <- data.frame(
        direction=rep(c("forward", "reverse"),3),
        proportion=c(forward_f1, reverse_f1,forward_f2, reverse_f2,forward_f3, reverse_f3),
        generation=c("F1", "F1", "F2", "F2", "F3", "F3")
        )


p <- ggplot(pdat, aes(x=direction, y=proportion)) +
  geom_col() +
  xlab("") +
  ylab("PO/GN gene number ratio") +
   facet_grid(rows=vars(generation)) +
   theme_bw() +
theme(
  strip.background = element_blank(),
  strip.text = element_blank(),
  legend.position = "none"
)

ggsave(file = "~/reciprocal_t/figures/forward_rev_ratios.pdf", p, w=2.5, h=6.3)
# 5.85 x 5.78


# get idea of how many of the total fall into the two categories

(length(forward_PO_f1) + length(forward_GN_f1) + length(reverse_PO_f1) + length(reverse_GN_f1))/length(f1_sig)

sum(f1_sig %in% unique(c(forward_PO_f1,forward_GN_f1,reverse_GN_f1,reverse_GN_f1)))/length(f1_sig)
sum(f2_sig %in% unique(c(forward_PO_f2,forward_GN_f2,reverse_GN_f2,reverse_GN_f2)))/length(f2_sig)
sum(f3_sig %in% unique(c(forward_PO_f3,forward_GN_f3,reverse_GN_f3,reverse_GN_f3)))/length(f3_sig)


######################
# plot points, across gens
######################

# make data frame for plotting
pdat <- data.frame(
        direction=rep(c("forward", "reverse"),6),
        count=c(length(forward_GN_f1), length(reverse_GN_f1), length(forward_PO_f1), length(reverse_PO_f1),
                     length(forward_GN_f2), length(reverse_GN_f2), length(forward_PO_f2), length(reverse_PO_f2),
                     length(forward_GN_f3), length(reverse_GN_f3), length(forward_PO_f3), length(reverse_PO_f3)),
        generation=c(rep("F1", 4), rep("F2", 4),rep("F3", 4)),
        class = rep(c("GN", "GN", "PO", "PO"),3)
        )
pdat

p <- ggplot(pdat, aes(x=direction, y=count, fill=class,color=class, shape=class, group=class)) +
  geom_line() +
  geom_point(size=4) +
  xlab("") +
  ylim(0,650)+
  ylab("Number of DE genes") +
   facet_grid(rows=vars(generation)) +
   theme_bw() +
theme(
  strip.background = element_blank(),
  strip.text = element_blank()
) +
  scale_fill_manual(values=c("grey80","black")) +
  scale_color_manual(values=c("grey80","black")) +
    scale_shape_manual(values = c(21,22))


ggsave(file = "~/reciprocal_t/figures/forward_rev_counts.pdf", p, w=2.5, h=3.5)




###########
# stats
###########
# following ho et al again.

library(DescTools)

f1 <- as.matrix(data.frame(
        PO = c(length(forward_PO_f1),length(reverse_PO_f1)),
        GN = c(length(forward_GN_f1),length(reverse_GN_f1))
        ))
row.names(f1) <- c("forward", "reverse")

GTest(f1) 

f2 <- as.matrix(data.frame(
        PO = c(length(forward_PO_f2),length(reverse_PO_f2)),
        GN = c(length(forward_GN_f2),length(reverse_GN_f2))
        ))
row.names(f2) <- c("forward", "reverse")

GTest(f2) 

f3 <- as.matrix(data.frame(
        PO = c(length(forward_PO_f3),length(reverse_PO_f3)),
        GN = c(length(forward_GN_f3),length(reverse_GN_f3))
        ))
row.names(f3) <- c("forward", "reverse")

GTest(f3)




#################################################
#################################################
#################################################
# venn
#################################################
#################################################
#################################################


df <- read.table("~/reciprocal_t/analysis/AF_change.txt", header=TRUE)

sum(df$sig_f1)
sum(df$sig_f2)
sum(df$sig_f3)

f1 <- df$SNP[which(df$sig_f1 == TRUE & df$sig_f2 == TRUE)
length(which(df$sig_f1 == TRUE & df$sig_f3 == TRUE))
length(which(df$sig_f2 == TRUE & df$sig_f3 == TRUE))
length(which(df$sig_f1 == TRUE & df$sig_f2 == TRUE& df$sig_f3 == TRUE))

# make venn

f1 <- as.character(df$SNP[which(df$sig_f1 == TRUE)])
f2 <- as.character(df$SNP[which(df$sig_f2 == TRUE)])
f3 <- as.character(df$SNP[which(df$sig_f3 == TRUE)])


library(VennDiagram)

v1 <- venn.diagram(x= list(F1=f1, F2=f2, F3=f3),
    filename=NULL,fill=c("dodgerblue3", "firebrick3", "springgreen3"),
        cat.cex = 1.5, cex=1.2)

pdf("~/reciprocal_t/figures/venn.pdf", h=5, w=5)
grid.newpage()
grid.draw(v1)
dev.off()

