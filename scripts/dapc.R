
library(stringr)
library(adegenet)
library(ggplot2)
library(ggpubr)
library(MCMCglmm)
library(dplyr)
library(tidyr)
library(DESeq2)
library(reshape)

########################################################################################################
########################################################################################################
########################################################################################################
##########################
########################## snps
##########################
########################################################################################################
########################################################################################################
########################################################################################################

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

dat <- af

###############################################################################################
###############################################################################################
######
###### F1
######
###############################################################################################
###############################################################################################
ctr <-dat[,grep("AAAA_F1|HHHH_F1",colnames(dat))] # in home envir
trt <-dat[,grep("AAHH_F1|HHAA_F1",colnames(dat))] # in away envir

pcp=prcomp(t(ctr), retx=TRUE, center=TRUE, scale.=FALSE)
scores=pcp$x
screeplot(pcp,bstick=T)

# adegenet: finding clusters (even though we know what clusters we want) -
clus.ctr=find.clusters(t(ctr),max.n.clus=10,
  n.pca=3, n.clust=2)

#Use clus$grp to rename groups
clus.ctr$grp=c(rep(2, 4),
                rep(1, 4))

# generate the actual DAPC
dp.ctr=dapc(t(ctr),clus.ctr$grp, var.loadings=TRUE, n.pca=3, n.da=1,
var.contrib =TRUE)
scatter(dp.ctr,bg="white",scree.da=FALSE,legend=TRUE,solid=.4)

# then add in the transplanted groups to the previously generated DAPC
pred.trt<-predict.dapc(dp.ctr,newdata=(t(trt)))

# create another dataframe structure in order to plot these predicted values
dpc=data.frame(rbind(dp.ctr$ind.coord,pred.trt$ind.scores))

dpc$Line <- substr(row.names(dpc),3,4)
dpc$Treatment <- substr(row.names(dpc),1,2)

dpc$gp <- c(rep(1.1,8), rep(1.9,4), rep(1.5,4))

# notusing this, but it can be a nice way to plot things in a density plot sort of way.
  # I like the other ones better though.
#z <- ggplot(dpc, aes(x=LD1, fill=Line, color=Treatment, linetype=Treatment)) +
#  geom_density(alpha=0.4, lwd=1) +
#  scale_color_manual(values=c("dodgerblue2","firebrick2"))+
#  scale_fill_manual(values=c("dodgerblue2","firebrick2")) +
#  geom_hline(color = "white", yintercept = 0, lwd=2.5) + theme_classic()

a <- ggplot(dpc, aes(x=LD1, y=gp, shape=Line, fill=Treatment)) +
  geom_point(size=6, lwd=2)+
  ylim(1, 3) +
  scale_color_manual(values=c("cornflowerblue","brown2"))+
  scale_fill_manual(values=alpha(c("cornflowerblue","brown2"), 0.8)) +
  scale_shape_manual(values=c(21, 24))+
  guides(fill=guide_legend(override.aes=list(shape=21, size=7, fill=c('cornflowerblue', 'brown2')),order = 2),
    shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)))) +
  geom_segment(
    x = mean(dpc$LD1[grep( "AAAA", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "HHAA", row.names(dpc))]),
    y=2.5, yend=2.5,
    size = 1.7, arrow = arrow(length = unit(0.2, "inches")),
    color="black")+
  geom_segment(
    x = mean(dpc$LD1[grep( "HHHH", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "AAHH", row.names(dpc))]),
    y=2.7, yend=2.7,
    size = 1.7, arrow = arrow(length = unit(0.2, "inches")),
    color="black")+
    theme_classic()+
  ylab(" ")+
  xlab(" ") +
  ggtitle("F1") +
  theme(plot.title = element_text(hjust = 0.5),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank())


##################
###### use MCMCglmm to set for sig
##################
dpc$rep <- substr(row.names(dpc), 9,13)

dpc$home <- 1
dpc$home[grep("AAAA|HHHH", row.names(dpc))] <- 0

prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V=1, nu=0.002,alpha.mu=0, alpha.V=1000)))

mod1 <- MCMCglmm(LD1~Line+Line:home, random=~rep, prior= prior, data=dpc,nitt=75000, thin=25, burnin=5000)
summary(mod1)
posterior.mode(mod1$VCV)

# calculating difference in magnitudes of orii:away and orio:away using sampled sets of parameters:
awayDelta=abs(mod1$Sol[,"LineHH:home"])-abs(mod1$Sol[,"LineAA:home"])

# 95% credible interval:
HPDinterval(awayDelta)
#         lower    upper
# var1 -1.870579 1.780305
if (is.na(table(awayDelta<0)[2])) {
  cat("p <",signif(1/length(awayDelta),1))
} else { cat("p =",signif(table(awayDelta<0)[2]/length(awayDelta),2)) } 
# p = 0.52

###############################################################################################
###############################################################################################
######
###### F2
######
###############################################################################################
###############################################################################################

ctr <-dat[,grep("AAAA_F2|HHHH_F2",colnames(dat))]
trt <-dat[,grep("AAHH_F2|HHAA_F2",colnames(dat))]

pcp=prcomp(t(ctr), retx=TRUE, center=TRUE, scale.=FALSE)
scores=pcp$x
screeplot(pcp,bstick=T)

# adegenet: finding clusters
clus.ctr=find.clusters(t(ctr),max.n.clus=10, 
  n.pca=3, n.clust=2)

#Use clus$grp to rename groups
clus.ctr$grp=c(rep(2, 4),
                rep(1, 4)) #tell the DF which groups you want to cluster; in this case in2in and off2off
# discriminant function for two groups:
dp.ctr=dapc(t(ctr),clus.ctr$grp, var.loadings=TRUE, n.pca=3, n.da=1,
var.contrib =TRUE) #keep 7 and 1
scatter(dp.ctr,bg="white",scree.da=FALSE,legend=TRUE,solid=.4) 

pred.trt<-predict.dapc(dp.ctr,newdata=(t(trt))) 

#must create another dataframe structure in order to plot these predicted values
dpc <- data.frame(rbind(dp.ctr$ind.coord,pred.trt$ind.scores))
dpc$Line <- substr(row.names(dpc),3,4)
dpc$Treatment <- substr(row.names(dpc),1,2)

dpc$gp <- c(rep(1.1,8), rep(1.9,4), rep(1.5,4))

b <- ggplot(dpc, aes(x=LD1, y=gp, shape=Line, fill=Treatment)) +
  geom_point(size=6, lwd=2)+ 
  ylim(1, 3) +
  scale_color_manual(values=c("cornflowerblue","brown2"))+
  scale_fill_manual(values=alpha(c("cornflowerblue","brown2"), 0.8)) +
  scale_shape_manual(values=c(21, 24))+
  guides(fill=guide_legend(override.aes=list(shape=21, size=7, fill=c('cornflowerblue', 'brown2')),order = 2),
    shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)))) +
  geom_segment(
    x = mean(dpc$LD1[grep( "AAAA", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "HHAA", row.names(dpc))]),
    y=2.5, yend=2.5,
     size = 1.7, arrow = arrow(length = unit(0.2, "inches")),
    color="black")+
  geom_segment(
    x = mean(dpc$LD1[grep( "HHHH", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "AAHH", row.names(dpc))]),
    y=2.7, yend=2.7,
     size = 1.7, arrow = arrow(length = unit(0.2, "inches")),
    color="black")+
    theme_classic()+
  ylab(" ")+
  xlab(" ") +
  ggtitle("F2") +
  theme(plot.title = element_text(hjust = 0.5),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank())

##################
###### use MCMCglmm to set for sig
##################
dpc$rep <- substr(row.names(dpc), 9,13)

dpc$home <- 1
dpc$home[grep("AAAA|HHHH", row.names(dpc))] <- 0

prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V=1, nu=0.002,alpha.mu=0, alpha.V=1000)))

mod1=MCMCglmm(LD1~Line+Line:home, random=~rep, prior= prior, data=dpc,nitt=75000, thin=25, burnin=5000)
#summary(mod1)
#posterior.mode(mod1$VCV)

# calculating difference in magnitudes of orii:away and orio:away using sampled sets of parameters:
awayDelta=abs(mod1$Sol[,"LineHH:home"])-abs(mod1$Sol[,"LineAA:home"])

# 95% credible interval:
HPDinterval(awayDelta)
#         lower   upper
# var1 -1.186448 1.85045
if (is.na(table(awayDelta<0)[2])) {
  cat("p <",signif(1/length(awayDelta),1))
} else { cat("p =",signif(table(awayDelta<0)[2]/length(awayDelta),2)) } 
# p = 0.3



###############################################################################################
###############################################################################################
######
###### SNPS F3
######
###############################################################################################
###############################################################################################

ctr <-dat[,grep("AAAA_F3|HHHH_F3",colnames(dat))] # in home envir
trt <-dat[,grep("AAHH_F3|HHAA_F3",colnames(dat))] # in away envir

pcp=prcomp(t(ctr), retx=TRUE, center=TRUE, scale.=FALSE) #scale.=TRUE
scores=pcp$x
screeplot(pcp,bstick=T)

clus.ctr=find.clusters(t(ctr),max.n.clus=10,
  n.pca=3, n.clust=2)

#Use clus$grp to rename groups
clus.ctr$grp=c(rep(2, 4),
                rep(1, 4))
# discriminant function for two groups:
dp.ctr=dapc(t(ctr),clus.ctr$grp, var.loadings=TRUE, n.pca=3, n.da=1,
var.contrib =TRUE)

pred.trt<-predict.dapc(dp.ctr,newdata=(t(trt)))

dpc=data.frame(rbind(dp.ctr$ind.coord,pred.trt$ind.scores))
dpc$Line <- substr(row.names(dpc),3,4)
dpc$Treatment <- substr(row.names(dpc),1,2)
dpc$gp <- c(rep(1.1,8), rep(1.9,4), rep(1.5,4))

c <- ggplot(dpc, aes(x=LD1, y=gp, shape=Line, fill=Treatment)) +
  geom_point(size=6, lwd=2)+ 
  ylim(1, 3) +
  scale_color_manual(values=c("cornflowerblue","brown2"))+
  scale_fill_manual(values=alpha(c("cornflowerblue","brown2"), 0.8)) +
  scale_shape_manual(values=c(21, 24))+
  guides(fill=guide_legend(override.aes=list(shape=21, size=7, fill=c('cornflowerblue', 'brown2')),order = 2),
    shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)))) +
  geom_segment(
    x = mean(dpc$LD1[grep( "AAAA", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "HHAA", row.names(dpc))]),
    y=2.5, yend=2.5,
     size = 1.7, arrow = arrow(length = unit(0.2, "inches")),
    color="black")+
  geom_segment(
    x = mean(dpc$LD1[grep( "HHHH", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "AAHH", row.names(dpc))]),
    y=2.7, yend=2.7,
     size = 1.7, arrow = arrow(length = unit(0.2, "inches")),
    color="black")+
    theme_classic()+
  ylab(" ")+
  xlab("Discriminant Function") +
  ggtitle("F3") +
  theme(plot.title = element_text(hjust = 0.5),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank())
#c


dpc$rep <- substr(row.names(dpc), 9,13)

dpc$home <- 1
dpc$home[grep("AAAA|HHHH", row.names(dpc))] <- 0

prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V=1, nu=0.002,alpha.mu=0, alpha.V=1000)))

mod1=MCMCglmm(LD1~Line+Line:home, random=~rep, prior= prior, data=dpc,nitt=75000, thin=25, burnin=5000)
summary(mod1)
posterior.mode(mod1$VCV)

# calculating difference in magnitudes of orii:away and orio:away using sampled sets of parameters:
awayDelta=abs(mod1$Sol[,"LineHH:home"])-abs(mod1$Sol[,"LineAA:home"])

# 95% credible interval:
HPDinterval(awayDelta)
#         lower   upper
# var1 0.009322394 2.77267
if (is.na(table(awayDelta<0)[2])) {
  cat("p <",signif(1/length(awayDelta),1))
} else { cat("p =",signif(table(awayDelta<0)[2]/length(awayDelta),2)) } 
# 0.018


png("~/reciprocal_t/figures/dapc_snps_gen.png", height = 6, width = 5, units="in", res=300)
ggarrange(a,b,c, ncol=1, nrow=3, common.legend=TRUE)
dev.off()

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
################
################ RNAseq
################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################

library(tximportData)
library(tximport)


###############################################################################################
###############################################################################################
######
###### DGE F1
######
###############################################################################################
###############################################################################################

###
### DEseq2-
###

dir <- "~/reciprocal_t/analysis/salmon"
#list.files(dir)

samples <- read.table("~/reciprocal_t/analysis/sample_id.txt", header=FALSE)

# now point to quant files
files <- file.path(dir, samples$V1, "quant.sf")
names(files) <- samples$V1
all(file.exists(files))

#subset files to only include F1
files <- files[grep("F1", files)]

gene_tran <- read.table("/data/atonsa/Atonsa_gen_trans_agp_gff/Atonsa_transcript_to_gene", header=FALSE)
tx2gene <- data.frame(transcript=gene_tran$V2, gene=gene_tran$V1)

# use tximport to read in files
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

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

# double check that the following agrees
rownames(sampleTable) <- colnames(txi$counts)

# setting this up, row info is each transcript/gene
  # column is phenotypic data

rownames(sampleTable) <- colnames(txi$counts)

# import to DESeq2. this usese counts from tximport. from salmon. gene level
dds <- DESeqDataSetFromTximport(txi,
                colData=sampleTable,
                design = ~ Line + Treatment + Line:Treatment)

# remove rows where count < 10 in more than 90% of samples
keep <- apply(counts(dds), 1, function(x) {ifelse(length(which(x > 10)) > 13, TRUE, FALSE)})
dds <- dds[keep,]
nrow(dds)
#[1] 23323

# then rolog transform
rld<-rlog(dds,blind=TRUE)
head(assay(rld))

dat=as.data.frame(assay(rld))
colnames(dat)<-colnames(dds)

####################
#discriminant function analysis
####################

ctr <-dat[,grep("AAAA|HHHH",colnames(dat))] # in home envir
trt <-dat[,grep("AAHH|HHAA",colnames(dat))] # in away envir

pcp=prcomp(t(ctr), retx=TRUE, center=TRUE, scale.=FALSE) # note. the scale command doesnt have much impact
scores=pcp$x
screeplot(pcp,bstick=T)

# adegenet: finding clusters (even though we know what clusters we want) - choose 7 PCs and 2 groups
clus.ctr=find.clusters(t(ctr),max.n.clus=7,
              n.pca=7, n.clust=2)  #keep 7 and 2

#rename groups
clus.ctr$grp=c(rep(2, 4),
                rep(1, 4))

# discriminant function for two groups:
dp.ctr=dapc(t(ctr),clus.ctr$grp,
              n.pca=4, n.da=2)  #keep 4 and 2
scatter(dp.ctr,bg="white",scree.da=FALSE,legend=TRUE,solid=.4) #discriminant function for ORIGIN type expression
pred.trt<-predict.dapc(dp.ctr,newdata=(t(trt))) #skip IO11C for host b/c outlier sample in WGCNA

#must create another dataframe structure in order to plot these predicted values
dpc=data.frame(rbind(dp.ctr$ind.coord,pred.trt$ind.scores))
dpc$Line <- substr(row.names(dpc),3,4)
dpc$Treatment <- substr(row.names(dpc),1,2)
dpc$gp <- c(rep(1.1,8), rep(1.9,4), rep(1.5,4))

a2 <- ggplot(dpc, aes(x=LD1, y=gp, fill=Treatment, shape=Line)) +
  geom_point(size=6, lwd=2)+ 
  ylim(1, 3) +
  scale_color_manual(values=c("cornflowerblue","brown2"))+
  scale_fill_manual(values=alpha(c("cornflowerblue","brown2"), 0.8)) +
  scale_shape_manual(values=c(21, 24))+
  guides(fill=guide_legend(override.aes=list(shape=21, size=7, fill=c('cornflowerblue', 'brown2')),order = 2),
    shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)))) +
  geom_segment(
    x = mean(dpc$LD1[grep( "AAAA", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "HHAA", row.names(dpc))]),
    y=2.5, yend=2.5,
     size = 1.7, arrow = arrow(length = unit(0.2, "inches")),
    color="black")+
  geom_segment(
    x = mean(dpc$LD1[grep( "HHHH", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "AAHH", row.names(dpc))]),
    y=2.7, yend=2.7,
     size = 1.7, arrow = arrow(length = unit(0.2, "inches")),
    color="black") +
  theme_classic()+
  ylab(" ")+
  xlab(" ") +
  ggtitle("F1") +
  theme(plot.title = element_text(hjust = 0.5),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank())


### stats:

dpc.all <- dpc
dpc.all$rep <- substr(row.names(dpc.all), 9,13)

dpc.all$home <- 1
dpc.all$home[grep("AAAA|HHHH", row.names(dpc.all))] <- 0

prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V=1, nu=0.002,alpha.mu=0, alpha.V=1000)))

mod1=MCMCglmm(LD1~Line+Line:home, random=~rep, prior= prior, data=dpc.all,nitt=75000, thin=25, burnin=5000)
summary(mod1)
posterior.mode(mod1$VCV)

# $Sol is the posterior distribution of the fixed effect
#head(mod1$Sol)

HPDinterval(mod1$Sol)

# calculating difference in magnitudes of orii:away and orio:away using sampled sets of parameters:
#awayDelta=abs(mod1$Sol[,"LineAA:home"])-abs(mod1$Sol[,"LineHH:home"])
awayDelta=abs(mod1$Sol[,"LineHH:home"])-abs(mod1$Sol[,"LineAA:home"])

# 95% credible interval:
HPDinterval(awayDelta)

if (is.na(table(awayDelta<0)[2])) {
  cat("p <",signif(1/length(awayDelta),1))
    } else { 
    cat("p =",signif(table(awayDelta<0)[2]/length(awayDelta),2)) 
  } # 0.97



###############################################################################################
###############################################################################################
######
###### DGE F2
######
###############################################################################################
###############################################################################################


#locate the directory containing the files. 
dir <- "~/reciprocal_t/analysis/salmon"
list.files(dir)

# make vector to point to quant files
  # read in table with sample ids
samples <- read.table("~/reciprocal_t/analysis/sample_id.txt", header=FALSE)

# now point to quant files
files <- file.path(dir, samples$V1, "quant.sf")
names(files) <- samples$V1
all(file.exists(files))

#subset files to only include F1 
files <- files[grep("F2", files)]

# associate transcripts with gene ids
# We first make a data.frame called tx2gene with two columns: 1) transcript ID and 2) gene ID. col order is important
gene_tran <- read.table("/data/atonsa/Atonsa_gen_trans_agp_gff/Atonsa_transcript_to_gene", header=FALSE)

tx2gene <- data.frame(transcript=gene_tran$V2, gene=gene_tran$V1)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)

head(txi$counts)

dir <- "~/reciprocal_t/analysis/salmon"
list.files(dir)

samples <- read.table("~/reciprocal_t/analysis/sample_id.txt", header=FALSE)

# now point to quant files
files <- file.path(dir, samples$V1, "quant.sf")
names(files) <- samples$V1
all(file.exists(files))

#subset files to only include F2
files <- files[grep("F2", files)]

gene_tran <- read.table("/data/atonsa/Atonsa_gen_trans_agp_gff/Atonsa_transcript_to_gene", header=FALSE)
tx2gene <- data.frame(transcript=gene_tran$V2, gene=gene_tran$V1)

f1samp <- as.data.frame(samples$V1[grep("F2", samples$V1)])
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

rownames(sampleTable) <- colnames(txi$counts)

# setting this up, row info is each transcript/gene
  # columb is phenotypic data
rownames(sampleTable) <- colnames(txi$counts)

dds <- DESeqDataSetFromTximport(txi,
                colData=sampleTable,
                design = ~ Line + Treatment + Line:Treatment)

# remove rows where count < 10 in more than 90% of samples
keep <- apply(counts(dds), 1, function(x) {ifelse(length(which(x > 10)) > 13, TRUE, FALSE)})
dds <- dds[keep,]
nrow(dds)
#[1] 24881

# then rolog transform
rld<-rlog(dds,blind=TRUE)
head(assay(rld))

dat=as.data.frame(assay(rld))
colnames(dat)<-colnames(dds)


####################
#discriminant function analysis
####################

ctr <-dat[,grep("AAAA|HHHH",colnames(dat))] # in home envir
trt <-dat[,grep("AAHH|HHAA",colnames(dat))] # in away envir

pcp=prcomp(t(ctr), retx=TRUE, center=TRUE, scale.=TRUE)
scores=pcp$x
screeplot(pcp,bstick=T)

# adegenet: finding clusters (even though we know what clusters we want) - choose 7 PCs and 2 groups
clus.ctr=find.clusters(t(ctr),max.n.clus=7,
              n.pca=7, n.clust=2)  #keep 7 and 2

#Use clus$grp to rename to in2in and off2off -
clus.ctr$grp=c(rep(2, 4),
                rep(1, 4)) #tell the DF which groups you want to cluster; in this case in2in and off2off

# discriminant function for two groups:
dp.ctr=dapc(t(ctr),clus.ctr$grp,
        n.pca=3, n.da=2)  #keep 7 and 2

scatter(dp.ctr,bg="white",scree.da=FALSE,legend=TRUE,solid=.4) #discriminant function for ORIGIN type expression
pred.trt<-predict.dapc(dp.ctr,newdata=(t(trt))) #skip IO11C for host b/c outlier sample in WGCNA

dpc <- data.frame(rbind(dp.ctr$ind.coord,pred.trt$ind.scores))
dpc$Line <- substr(row.names(dpc),3,4)
dpc$Treatment <- substr(row.names(dpc),1,2)
dpc$gp <- c(rep(1.1,8), rep(1.9,4), rep(1.5,4))

b2 <- ggplot(dpc, aes(x=LD1, y=gp, fill=Treatment, shape=Line)) +
  geom_point(size=6, lwd=2)+ 
  ylim(1, 3) +
  scale_color_manual(values=c("cornflowerblue","brown2"))+
  scale_fill_manual(values=alpha(c("cornflowerblue","brown2"), 0.8)) +
  scale_shape_manual(values=c(21, 24))+
  guides(fill=guide_legend(override.aes=list(shape=21, size=7, fill=c('cornflowerblue', 'brown2')),order = 2),
    shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)))) +
  geom_segment(
    x = mean(dpc$LD1[grep( "AAAA", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "HHAA", row.names(dpc))]),
    y=2.5, yend=2.5,
     size = 1.7, arrow = arrow(length = unit(0.2, "inches")),
    color="black")+
  geom_segment(
    x = mean(dpc$LD1[grep( "HHHH", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "AAHH", row.names(dpc))]),
    y=2.7, yend=2.7,
     size = 1.7, arrow = arrow(length = unit(0.2, "inches")),
    color="black")+
  theme_classic()+
  ylab(" ")+
  xlab(" ") +
  ggtitle("F2") +
  theme(plot.title = element_text(hjust = 0.5),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank())

### stats:
dpc.all <- dpc
dpc.all$rep <- substr(row.names(dpc.all), 9,13)

dpc.all$home <- 1
dpc.all$home[grep("AAAA|HHHH", row.names(dpc.all))] <- 0

prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V=1, nu=0.002,alpha.mu=0, alpha.V=1000)))

mod1=MCMCglmm(LD1~Line+Line:home, random=~rep, prior= prior, data=dpc.all,nitt=75000, thin=25, burnin=5000)
summary(mod1)
posterior.mode(mod1$VCV)

# $Sol is the posterior distribution of the fixed effect
#head(mod1$Sol)
HPDinterval(mod1$Sol)

# calculating difference in magnitudes of orii:away and orio:away using sampled sets of parameters:
#awayDelta=abs(mod1$Sol[,"LineAA:home"])-abs(mod1$Sol[,"LineHH:home"])
awayDelta=abs(mod1$Sol[,"LineHH:home"])-abs(mod1$Sol[,"LineAA:home"])

# 95% credible interval:
HPDinterval(awayDelta)

if (is.na(table(awayDelta<0)[2])) {
  cat("p <",signif(1/length(awayDelta),1))
    } else {
    cat("p =",signif(table(awayDelta<0)[2]/length(awayDelta),2))
  } # 0.81

###############################################################################################
###############################################################################################
######
###### DGE F3
######
###############################################################################################
###############################################################################################

samples <- read.table("~/reciprocal_t/analysis/sample_id.txt", header=FALSE)

# now point to quant files
files <- file.path(dir, samples$V1, "quant.sf")
names(files) <- samples$V1
all(file.exists(files))

#subset files to only include F3
files <- files[grep("F3", files)]

gene_tran <- read.table("/data/atonsa/Atonsa_gen_trans_agp_gff/Atonsa_transcript_to_gene", header=FALSE)

tx2gene <- data.frame(transcript=gene_tran$V2, gene=gene_tran$V1)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)

head(txi$counts)

tx2gene <- data.frame(transcript=gene_tran$V2, gene=gene_tran$V1)

f1samp <- as.data.frame(samples$V1[grep("F3", samples$V1)])
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

# double check that the following agrees
rownames(sampleTable) <- colnames(txi$counts)

dds <- DESeqDataSetFromTximport(txi, 
                colData=sampleTable, 
                design = ~ Line + Treatment + Line:Treatment)

# remove rows where count < 10 in more than 90% of samples
keep <- apply(counts(dds), 1, function(x) {ifelse(length(which(x > 10)) > 13, TRUE, FALSE)})
dds <- dds[keep,]
nrow(dds)
#[1] 24131

# then rolog transform
rld<-rlog(dds,blind=TRUE)

dat=as.data.frame(assay(rld))
colnames(dat)<-colnames(dds)

####################
#discriminant function analysis
####################

ctr <-dat[,grep("AAAA|HHHH",colnames(dat))] # in home envir
trt <-dat[,grep("AAHH|HHAA",colnames(dat))] # in away envir

pcp=prcomp(t(ctr), retx=TRUE, center=TRUE, scale.=TRUE)
scores=pcp$x
screeplot(pcp,bstick=T)

# adegenet: finding clusters (even though we know what clusters we want) - choose 7 PCs and 2 groups
clus.ctr=find.clusters(t(ctr),max.n.clus=7,
             n.pca=7, n.clust=2)  #keep 7 and 2


#Use clus$grp to rename to in2in and off2off -
clus.ctr$grp=c(rep(2, 4),
                rep(1, 4)) #tell the DF which groups you want to cluster; in this case in2in and off2off

# discriminant function for two groups:
dp.ctr=dapc(t(ctr),clus.ctr$grp,
             n.pca=4, n.da=2)  #keep 7 and 2

scatter(dp.ctr,bg="white",scree.da=FALSE,legend=TRUE,solid=.4) #discriminant function for ORIGIN type expression
pred.trt<-predict.dapc(dp.ctr,newdata=(t(trt))) #skip IO11C for host b/c outlier sample in WGCNA

#assign groups
dpc=data.frame(rbind(dp.ctr$ind.coord,pred.trt$ind.scores))
dpc$Line <- substr(row.names(dpc),3,4)
dpc$Treatment <- substr(row.names(dpc),1,2)
dpc$gp <- c(rep(1.1,8), rep(1.9,4), rep(1.5,4))

c2 <- ggplot(dpc, aes(x=LD1, y=gp, fill=Treatment, shape=Line)) +
  geom_point(size=6, lwd=2)+ 
  ylim(1, 3) +
  scale_color_manual(values=c("cornflowerblue","brown2"))+
  scale_fill_manual(values=alpha(c("cornflowerblue","brown2"), 0.8)) +
  scale_shape_manual(values=c(21, 24))+
  guides(fill=guide_legend(override.aes=list(shape=21, size=7, fill=c('cornflowerblue', 'brown2')),order = 2),
    shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)))) +
  geom_segment(
    x = mean(dpc$LD1[grep( "AAAA", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "HHAA", row.names(dpc))]),
    y=2.5, yend=2.5,
     size = 1.7, arrow = arrow(length = unit(0.2, "inches")),
    color="black")+
  geom_segment(
    x = mean(dpc$LD1[grep( "HHHH", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "AAHH", row.names(dpc))]),
    y=2.7, yend=2.7,
     size = 1.7, arrow = arrow(length = unit(0.2, "inches")),
    color="black")+
  theme_classic()+
  ylab(" ")+
  xlab("Discriminant Function") +
  ggtitle("F3") +
  theme(plot.title = element_text(hjust = 0.5),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank())

png("~/reciprocal_t/figures/dapc_RNA_gen.png", height = 6, width = 5, units="in", res=300)
ggarrange(a2,b2,c2, ncol=1, nrow=3, common.legend=TRUE)
dev.off()

### stats:

dpc.all <- dpc
dpc.all$rep <- substr(row.names(dpc.all), 9,13)

dpc.all$home <- 1
dpc.all$home[grep("AAAA|HHHH", row.names(dpc.all))] <- 0

prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V=1, nu=0.002,alpha.mu=0, alpha.V=1000)))

mod1=MCMCglmm(LD1~Line+Line:home, random=~rep, prior= prior, data=dpc.all,nitt=75000, thin=25, burnin=5000)
summary(mod1)
posterior.mode(mod1$VCV)

HPDinterval(mod1$Sol)

# calculating difference in magnitudes of orii:away and orio:away using sampled sets of parameters:
#awayDelta=abs(mod1$Sol[,"LineAA:home"])-abs(mod1$Sol[,"LineHH:home"])
awayDelta=abs(mod1$Sol[,"LineHH:home"])-abs(mod1$Sol[,"LineAA:home"])

# 95% credible interval:
HPDinterval(awayDelta)

if (is.na(table(awayDelta<0)[2])) {
  cat("p <",signif(1/length(awayDelta),1))
    } else {
    cat("p =",signif(table(awayDelta<0)[2]/length(awayDelta),2))
  } #0.61


library("ggpubr")

png("~/reciprocal_t/figures/dapc_all_gen1.png", height = 6, width = 8, units="in", res=300)
ggarrange(a, a2, b, b2, c, c2, ncol=2, nrow=3, common.legend=TRUE)
dev.off()
