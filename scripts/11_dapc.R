
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

gene_tran <- read.table("/data/copepods/tonsa_transcriptome/Atonsa_gen_trans_agp_gff/Atonsa_transcript_to_gene", header=FALSE)
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
#head(assay(rld))

dat=as.data.frame(assay(rld))
colnames(dat)<-colnames(dds)

####################
#discriminant function analysis
####################

ctr <-dat[,grep("AAAA|HHHH",colnames(dat))] # in home envir
trt <-dat[,grep("AAHH|HHAA",colnames(dat))] # in away envir

pcp=prcomp(t(ctr), retx=TRUE, center=TRUE, scale.=FALSE) # note. the scale command doesnt have much impact
scores=pcp$x
#screeplot(pcp,bstick=T)

# adegenet: finding clusters (even though we know what clusters we want) - choose 7 PCs and 2 groups
clus.ctr=find.clusters(t(ctr),max.n.clus=7,
              n.pca=7, n.clust=2)  #keep 7 and 2

#rename groups
clus.ctr$grp=c(rep(2, 4),
                rep(1, 4))

# discriminant function for two groups:
dp.ctr=dapc(t(ctr),clus.ctr$grp,
              n.pca=4, n.da=2)  #keep 4 and 2
#scatter(dp.ctr,bg="white",scree.da=FALSE,legend=TRUE,solid=.4) #discriminant function for ORIGIN type expression
pred.trt<-predict.dapc(dp.ctr,newdata=(t(trt))) #skip IO11C for host b/c outlier sample in WGCNA

#must create another dataframe structure in order to plot these predicted values
dpc=data.frame(rbind(dp.ctr$ind.coord,pred.trt$ind.scores))
dpc$Line <- substr(row.names(dpc),3,4)
dpc$Treatment <- substr(row.names(dpc),1,2)
dpc$gp <- c(rep(1.1,4), rep(1.5,8),rep(1.1,4))
dpc$LD1 <- dpc$LD1*-1

dpc.dge.1 <- dpc

gp_means <- dpc %>% group_by(Line, Treatment) %>% 
 summarize(mean=mean(LD1))
gp_means$gp <- c(1.1, 1.1,1.5, 1.5)

a2 <- ggplot(dpc, aes(x=LD1, y=gp, fill=Treatment, shape=Line)) +
  geom_point(size=3, lwd=2)+ 
  ylim(1, 2) +
  xlim(-5,5) +
  scale_color_manual(values=c("#6699CC","#CC3333"))+
  scale_fill_manual(values=c("#6699CC","#CC3333")) +
  scale_shape_manual(values=c(21, 24))+
  guides(fill=guide_legend(override.aes=list(shape=21, size=7, fill=c("#6699CC","#CC3333")),order = 2),
    shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)))) +
  #stat_summary(fun="mean") +
  theme_classic()+
  ylab(" ")+
  xlab(" ") +
  ggtitle("F1") +
  theme(plot.title = element_text(hjust = 0.5),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank()) +
  geom_point(data=gp_means, aes(x=mean, y=gp, fill=Treatment, shape=Line), size=7, alpha=0.5) +
  geom_segment(
    x = mean(dpc$LD1[grep( "AAAA", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "HHAA", row.names(dpc))]),
    y=1.1, yend=1.1,
    lty=1,
     size = 1, arrow = arrow(length = unit(0.2, "inches")),
    color="black")+
  geom_segment(
    x = mean(dpc$LD1[grep( "HHHH", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "AAHH", row.names(dpc))]),
    y=1.5, yend=1.5,
     size = 1, arrow = arrow(length = unit(0.2, "inches")),
    color="gray50")+
  theme(legend.position = "none")


### stats:

dpc.all <- dpc
dpc.all$rep <- substr(row.names(dpc.all), 9,13)

dpc.all$home <- 1
dpc.all$home[grep("AAAA|HHHH", row.names(dpc.all))] <- 0

prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V=1, nu=0.002,alpha.mu=0, alpha.V=1000)))

mod1=MCMCglmm(LD1~Line+Line:home, random=~rep, prior= prior, data=dpc.all,nitt=75000, thin=25, burnin=5000)
summary(mod1)
#            post.mean l-95% CI u-95% CI eff.samp   pMCMC
#(Intercept)   -3.2919  -4.5618  -2.0473     2800 0.00286 **
#LineHH         6.6354   5.5284   7.7635     3059 < 4e-04 ***
#LineAA:home    3.5265   2.4352   4.6981     2800 < 4e-04 ***
#LineHH:home   -1.9568  -3.0912  -0.9127     2836 0.00143 **
posterior.mode(mod1$VCV)

summary(mod1$Sol)
# $Sol is the posterior distribution of the fixed effect
#head(mod1$Sol)

HPDinterval(mod1$Sol)

# calculating difference in magnitudes of LineAA:home and LineHH:home using sampled sets of parameters:
awayDelta=abs(mod1$Sol[,"LineAA:home"]) -abs(mod1$Sol[,"LineHH:home"])

# 95% credible interval:
HPDinterval(awayDelta)

if (is.na(table(awayDelta<0)[2])) {
  cat("p <",signif(1/length(awayDelta),1))
    } else { 
    cat("p =",signif(table(awayDelta<0)[2]/length(awayDelta),2)) 
  }


###############################################################################################
###############################################################################################
######
###### DGE F2
######
###############################################################################################
###############################################################################################

# now point to quant files
files <- file.path(dir, samples$V1, "quant.sf")
names(files) <- samples$V1
all(file.exists(files))

#subset files to only include F1 
files <- files[grep("F2", files)]

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)


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


dat=as.data.frame(assay(rld))
colnames(dat)<-colnames(dds)


####################
#discriminant function analysis
####################

ctr <-dat[,grep("AAAA|HHHH",colnames(dat))] # in home envir
trt <-dat[,grep("AAHH|HHAA",colnames(dat))] # in away envir

pcp=prcomp(t(ctr), retx=TRUE, center=TRUE, scale.=TRUE)
scores=pcp$x
#screeplot(pcp,bstick=T)

# adegenet: finding clusters (even though we know what clusters we want) - choose 7 PCs and 2 groups
clus.ctr=find.clusters(t(ctr),max.n.clus=7,
              n.pca=7, n.clust=2)  #keep 7 and 2

#Use clus$grp to rename to in2in and off2off -
clus.ctr$grp=c(rep(2, 4),
                rep(1, 4)) #tell the DF which groups you want to cluster; in this case in2in and off2off

# discriminant function for two groups:
dp.ctr=dapc(t(ctr),clus.ctr$grp,
        n.pca=3, n.da=2)  #keep 7 and 2

#scatter(dp.ctr,bg="white",scree.da=FALSE,legend=TRUE,solid=.4) #discriminant function for ORIGIN type expression
pred.trt<-predict.dapc(dp.ctr,newdata=(t(trt))) #skip IO11C for host b/c outlier sample in WGCNA

dpc <- data.frame(rbind(dp.ctr$ind.coord,pred.trt$ind.scores))
dpc$Line <- substr(row.names(dpc),3,4)
dpc$Treatment <- substr(row.names(dpc),1,2)
dpc$gp <- c(rep(1.1,4), rep(1.5,8),rep(1.1,4))
dpc$LD1 <- dpc$LD1*-1
dpc.dge.2 <- dpc


gp_means <- dpc %>% group_by(Line, Treatment) %>% 
 summarize(mean=mean(LD1))
gp_means$gp <- c(1.1, 1.1,1.5, 1.5)

b2 <- ggplot(dpc, aes(x=LD1, y=gp, fill=Treatment, shape=Line)) +
  geom_point(size=3, lwd=2)+ 
  ylim(1, 2) +
  xlim(-5, 5) +
  scale_color_manual(values=c("#6699CC","#CC3333"))+
  scale_fill_manual(values=c("#6699CC","#CC3333")) +
  scale_shape_manual(values=c(21, 24))+
  guides(fill=guide_legend(override.aes=list(shape=21, size=7, fill=c("#6699CC","#CC3333")),order = 2),
    shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)))) +
  #stat_summary(fun="mean") +
  theme_classic()+
  ylab(" ")+
  xlab(" ") +
  ggtitle("F2") +
  theme(plot.title = element_text(hjust = 0.5),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank()) +
  geom_point(data=gp_means, aes(x=mean, y=gp, fill=Treatment, shape=Line), size=7, alpha=0.5) +
  geom_segment(
    x = mean(dpc$LD1[grep( "AAAA", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "HHAA", row.names(dpc))]),
    y=1.1, yend=1.1,
    lty=1,
     size = 1, arrow = arrow(length = unit(0.2, "inches")),
    color="gray50")+
  geom_segment(
    x = mean(dpc$LD1[grep( "HHHH", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "AAHH", row.names(dpc))]),
    y=1.5, yend=1.5,
     size = 1, arrow = arrow(length = unit(0.2, "inches")),
    color="gray50")+
  theme(legend.position = "none")

### stats:
dpc.all <- dpc
dpc.all$rep <- substr(row.names(dpc.all), 9,13)

dpc.all$home <- 1
dpc.all$home[grep("AAAA|HHHH", row.names(dpc.all))] <- 0

prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V=1, nu=0.002,alpha.mu=0, alpha.V=1000)))

mod1=MCMCglmm(LD1~Line+Line:home, random=~rep, prior= prior, data=dpc.all,nitt=75000, thin=25, burnin=5000)
summary(mod1)
#             post.mean l-95% CI u-95% CI eff.samp    pMCMC
#(Intercept)    3.3268   2.1485   4.5441     2652 0.001429 **
#LineHH        -6.6198  -7.7842  -5.4205     3071  < 4e-04 ***
#LineAA:home   -3.5292  -4.6578  -2.4419     2800 0.000714 ***
#LineHH:home    1.9427   0.8655   3.0950     2800 0.003571 **
posterior.mode(mod1$VCV)

# $Sol is the posterior distribution of the fixed effect
#head(mod1$Sol)
HPDinterval(mod1$Sol)

awayDelta=abs(mod1$Sol[,"LineAA:home"])-abs(mod1$Sol[,"LineHH:home"])

# 95% credible interval:
HPDinterval(awayDelta)

if (is.na(table(awayDelta<0)[2])) {
  cat("p <",signif(1/length(awayDelta),1))
    } else {
    cat("p =",signif(table(awayDelta<0)[2]/length(awayDelta),2))
  }

###############################################################################################
###############################################################################################
######
###### DGE F3
######
###############################################################################################
###############################################################################################

# now point to quant files
files <- file.path(dir, samples$V1, "quant.sf")
names(files) <- samples$V1
all(file.exists(files))

#subset files to only include F3
files <- files[grep("F3", files)]

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

# adegenet: finding clusters (even though we know what clusters we want) - choose 7 PCs and 2 groups
clus.ctr=find.clusters(t(ctr),max.n.clus=7,
             n.pca=7, n.clust=2)  #keep 7 and 2


#Use clus$grp to rename to in2in and off2off -
clus.ctr$grp=c(rep(2, 4),
                rep(1, 4)) #tell the DF which groups you want to cluster; in this case in2in and off2off

# discriminant function for two groups:
dp.ctr=dapc(t(ctr),clus.ctr$grp,
             n.pca=4, n.da=2)  #keep 7 and 2

#scatter(dp.ctr,bg="white",scree.da=FALSE,legend=TRUE,solid=.4) #discriminant function for ORIGIN type expression
pred.trt<-predict.dapc(dp.ctr,newdata=(t(trt))) #skip IO11C for host b/c outlier sample in WGCNA

#assign groups
dpc <- data.frame(rbind(dp.ctr$ind.coord,pred.trt$ind.scores))
dpc$Line <- substr(row.names(dpc),3,4)
dpc$Treatment <- substr(row.names(dpc),1,2)
dpc$gp <- c(rep(1.1,4), rep(1.5,8),rep(1.1,4))
dpc$LD1 <- dpc$LD1*-1

dpc.dge.3 <- dpc


dpc_dge <- rbind(dpc.dge.1, dpc.dge.2, dpc.dge.3)

write.table(dpc_dge, "~/reciprocal_t/analysis/dpc_dge.txt", sep="\t", quote=F, row.names=F)






gp_means <- dpc %>% group_by(Line, Treatment) %>% 
 summarize(mean=mean(LD1))
gp_means$gp <- c(1.1, 1.1,1.5, 1.5)

c2 <- ggplot(dpc, aes(x=LD1, y=gp, fill=Treatment, shape=Line)) +
  geom_point(size=3, lwd=2)+ 
  ylim(1, 2) +
  xlim(-5, 5) +
  scale_color_manual(values=c("#6699CC","#CC3333"))+
  scale_fill_manual(values=c("#6699CC","#CC3333")) +
  scale_shape_manual(values=c(21, 24))+
  guides(fill=guide_legend(override.aes=list(shape=21, size=7, fill=c("#6699CC","#CC3333")),order = 2),
    shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)))) +
  #stat_summary(fun="mean") +
  theme_classic()+
  ylab(" ")+
  xlab(" ") +
  ggtitle("F3") +
  theme(plot.title = element_text(hjust = 0.5),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank()) +
  geom_point(data=gp_means, aes(x=mean, y=gp, fill=Treatment, shape=Line), size=7, alpha=0.5) +
  geom_segment(
    x = mean(dpc$LD1[grep( "AAAA", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "HHAA", row.names(dpc))]),
    y=1.1, yend=1.1,
    lty=1,
     size = 1, arrow = arrow(length = unit(0.2, "inches")),
    color="gray50")+
  geom_segment(
    x = mean(dpc$LD1[grep( "HHHH", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "AAHH", row.names(dpc))]),
    y=1.5, yend=1.5,
     size = 1, arrow = arrow(length = unit(0.2, "inches")),
    color="gray50")+
  theme(legend.position = "none")

pdf("~/Documents/UVM/Reciprocal_transplant/figures/dapc_RNA_gen.pdf", height = 6, width = 4)
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
#            post.mean l-95% CI u-95% CI eff.samp  pMCMC
#(Intercept)    -3.480   -4.805   -2.140     2800 <4e-04 ***
#LineHH          6.951    5.492    8.494     2800 <4e-04 ***
#LineAA:home     4.799    3.306    6.291     2800 <4e-04 ***
#LineHH:home    -4.555   -6.071   -3.047     2800 <4e-04 ***
posterior.mode(mod1$VCV)

HPDinterval(mod1$Sol)

#awayDelta=abs(mod1$Sol[,"LineAA:home"])-abs(mod1$Sol[,"LineHH:home"])
awayDelta=abs(mod1$Sol[,"LineHH:home"])-abs(mod1$Sol[,"LineAA:home"])

# 95% credible interval:
HPDinterval(awayDelta)

if (is.na(table(awayDelta<0)[2])) {
  cat("p <",signif(1/length(awayDelta),1))
    } else {
    cat("p =",signif(table(awayDelta<0)[2]/length(awayDelta),2))
  } 






# all gens:
#f1
#LineAA:home    3.5265   2.4352   4.6981     2800 < 4e-04 ***
#LineHH:home   -1.9568  -3.0912  -0.9127     2836 0.00143 **
#f2
#LineAA:home   -3.5292  -4.6578  -2.4419     2800 0.000714 ***
#LineHH:home    1.9427   0.8655   3.0950     2800 0.003571 **
# f3
#LineAA:home     4.799    3.306    6.291     2800 <4e-04 ***
#LineHH:home    -4.555   -6.071   -3.047     2800 <4e-04 ***

p.adjust(c(4e-04,0.00143,0.000714,0.003571,4e-04,4e-04 ), method="bonferroni")
#[1] 0.002400 0.008580 0.004284 0.021426 0.002400 0.002400



######################################################################################################################
###########################################################
###########################################################
## allele freqs
###########################################################
###########################################################
######################################################################################################################


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
#screeplot(pcp,bstick=T)

# adegenet: finding clusters (even though we know what clusters we want) -
clus.ctr=find.clusters(t(ctr),max.n.clus=10,
  n.pca=3, n.clust=2)

#Use clus$grp to rename groups
clus.ctr$grp=c(rep(2, 4),
                rep(1, 4))

# generate the actual DAPC
dp.ctr=dapc(t(ctr),clus.ctr$grp, var.loadings=TRUE, n.pca=3, n.da=1,
var.contrib =TRUE)
#scatter(dp.ctr,bg="white",scree.da=FALSE,legend=TRUE,solid=.4)

# then add in the transplanted groups to the previously generated DAPC
pred.trt<-predict.dapc(dp.ctr,newdata=(t(trt)))

# create another dataframe structure in order to plot these predicted values

dpc=data.frame(rbind(dp.ctr$ind.coord,pred.trt$ind.scores))
dpc$Line <- substr(row.names(dpc),3,4)
dpc$Treatment <- substr(row.names(dpc),1,2)
dpc$gp <- c(rep(1.1,4), rep(1.5,8),rep(1.1,4))
dpc$LD1 <- dpc$LD1*-1

dpc1 <- dpc
gp_means <- dpc %>% group_by(Line, Treatment) %>% 
 summarize(mean=mean(LD1))
gp_means$gp <- c(1.1, 1.1,1.5, 1.5)


a2 <- ggplot(dpc, aes(x=LD1, y=gp, fill=Treatment, shape=Line)) +
  geom_point(size=3, lwd=2)+ 
  ylim(1, 2) +
  xlim(-5.5,5.5) +
  scale_color_manual(values=c("#6699CC","#CC3333"))+
  scale_fill_manual(values=c("#6699CC","#CC3333")) +
  scale_shape_manual(values=c(21, 24))+
  guides(fill=guide_legend(override.aes=list(shape=21, size=7, fill=c("#6699CC","#CC3333")),order = 2),
    shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)))) +
  #stat_summary(fun="mean") +
  theme_classic()+
  ylab(" ")+
  xlab(" ") +
  ggtitle("F1") +
  theme(plot.title = element_text(hjust = 0.5),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank()) +
  geom_point(data=gp_means, aes(x=mean, y=gp, fill=Treatment, shape=Line), size=7, alpha=0.5) +
  geom_segment(
    x = mean(dpc$LD1[grep( "AAAA", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "HHAA", row.names(dpc))]),
    y=1.1, yend=1.1,
    lty=1,
     size = 1, arrow = arrow(length = unit(0.2, "inches")),
    color="gray50")+
  geom_segment(
    x = mean(dpc$LD1[grep( "HHHH", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "AAHH", row.names(dpc))]),
    y=1.5, yend=1.5,
     size = 1, arrow = arrow(length = unit(0.2, "inches")),
    color="gray50")+
  theme(legend.position = "none")

###############################################################################################
###############################################################################################
######
###### F2
######
###############################################################################################
###############################################################################################
ctr <-dat[,grep("AAAA_F2|HHHH_F2",colnames(dat))] # in home envir
trt <-dat[,grep("AAHH_F2|HHAA_F2",colnames(dat))] # in away envir

pcp=prcomp(t(ctr), retx=TRUE, center=TRUE, scale.=FALSE)
scores=pcp$x
#screeplot(pcp,bstick=T)

# adegenet: finding clusters (even though we know what clusters we want) -
clus.ctr=find.clusters(t(ctr),max.n.clus=10,
  n.pca=3, n.clust=2)

#Use clus$grp to rename groups
clus.ctr$grp=c(rep(2, 4),
                rep(1, 4))

# generate the actual DAPC
dp.ctr=dapc(t(ctr),clus.ctr$grp, var.loadings=TRUE, n.pca=3, n.da=1,
var.contrib =TRUE)
#scatter(dp.ctr,bg="white",scree.da=FALSE,legend=TRUE,solid=.4)

# then add in the transplanted groups to the previously generated DAPC
pred.trt<-predict.dapc(dp.ctr,newdata=(t(trt)))

# create another dataframe structure in order to plot these predicted values
dpc=data.frame(rbind(dp.ctr$ind.coord,pred.trt$ind.scores))
dpc$Line <- substr(row.names(dpc),3,4)
dpc$Treatment <- substr(row.names(dpc),1,2)
dpc$gp <- c(rep(1.1,4), rep(1.5,8),rep(1.1,4))
dpc$LD1 <- dpc$LD1*-1
dpc2 <- dpc

gp_means <- dpc %>% group_by(Line, Treatment) %>% 
 summarize(mean=mean(LD1))
gp_means$gp <- c(1.1, 1.1,1.5, 1.5)

b2 <- ggplot(dpc, aes(x=LD1, y=gp, fill=Treatment, shape=Line)) +
  geom_point(size=3, lwd=2)+ 
  ylim(1, 2) +
  xlim(-5.5,5.5) +
  scale_color_manual(values=c("#6699CC","#CC3333"))+
  scale_fill_manual(values=c("#6699CC","#CC3333")) +
  scale_shape_manual(values=c(21, 24))+
  guides(fill=guide_legend(override.aes=list(shape=21, size=7, fill=c("#6699CC","#CC3333")),order = 2),
    shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)))) +
  #stat_summary(fun="mean") +
  theme_classic()+
  ylab(" ")+
  xlab(" ") +
  ggtitle("F2") +
  theme(plot.title = element_text(hjust = 0.5),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank()) +
  geom_point(data=gp_means, aes(x=mean, y=gp, fill=Treatment, shape=Line), size=7, alpha=0.5) +
  geom_segment(
    x = mean(dpc$LD1[grep( "AAAA", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "HHAA", row.names(dpc))]),
    y=1.1, yend=1.1,
    lty=1,
     size = 1, arrow = arrow(length = unit(0.2, "inches")),
    color="gray50")+
  geom_segment(
    x = mean(dpc$LD1[grep( "HHHH", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "AAHH", row.names(dpc))]),
    y=1.5, yend=1.5,
     size = 1, arrow = arrow(length = unit(0.2, "inches")),
    color="gray50")+
  theme(legend.position = "none")

  ###############################################################################################
###############################################################################################
######
###### F3
######
###############################################################################################
###############################################################################################
ctr <-dat[,grep("AAAA_F3|HHHH_F3",colnames(dat))] # in home envir
trt <-dat[,grep("AAHH_F3|HHAA_F3",colnames(dat))] # in away envir

pcp=prcomp(t(ctr), retx=TRUE, center=TRUE, scale.=FALSE)
scores=pcp$x
#screeplot(pcp,bstick=T)

# adegenet: finding clusters (even though we know what clusters we want) -
clus.ctr=find.clusters(t(ctr),max.n.clus=10,
  n.pca=3, n.clust=2)

#Use clus$grp to rename groups
clus.ctr$grp=c(rep(2, 4),
                rep(1, 4))

# generate the actual DAPC
dp.ctr=dapc(t(ctr),clus.ctr$grp, var.loadings=TRUE, n.pca=3, n.da=1,
var.contrib =TRUE)
#scatter(dp.ctr,bg="white",scree.da=FALSE,legend=TRUE,solid=.4)

# then add in the transplanted groups to the previously generated DAPC
pred.trt<-predict.dapc(dp.ctr,newdata=(t(trt)))

# create another dataframe structure in order to plot these predicted values
dpc=data.frame(rbind(dp.ctr$ind.coord,pred.trt$ind.scores))
dpc=data.frame(rbind(dp.ctr$ind.coord,pred.trt$ind.scores))
dpc$Line <- substr(row.names(dpc),3,4)
dpc$Treatment <- substr(row.names(dpc),1,2)
dpc$gp <- c(rep(1.1,4), rep(1.5,8),rep(1.1,4))
dpc$LD1 <- dpc$LD1*-1
dpc3 <- dpc


dpc_af <- rbind(dpc1, dpc2, dpc3)

write.table(dpc_af, "~/reciprocal_t/analysis/dpc_af.txt", sep="\t", quote=F, row.names=F)

gp_means <- dpc %>% group_by(Line, Treatment) %>% 
 summarize(mean=mean(LD1))
gp_means$gp <- c(1.1, 1.1,1.5, 1.5)

c2 <- ggplot(dpc, aes(x=LD1, y=gp, fill=Treatment, shape=Line)) +
  geom_point(size=3, lwd=2)+ 
  ylim(1, 2) +
  xlim(-5.5,5.5) +
  scale_color_manual(values=c("#6699CC","#CC3333"))+
  scale_fill_manual(values=c("#6699CC","#CC3333")) +
  scale_shape_manual(values=c(21, 24))+
  guides(fill=guide_legend(override.aes=list(shape=21, size=7, fill=c("#6699CC","#CC3333")),order = 2),
    shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)))) +
  #stat_summary(fun="mean") +
  theme_classic()+
  ylab(" ")+
  xlab(" ") +
  ggtitle("F3") +
  theme(plot.title = element_text(hjust = 0.5),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank()) +
  geom_point(data=gp_means, aes(x=mean, y=gp, fill=Treatment, shape=Line), size=7, alpha=0.5) +
  geom_segment(
    x = mean(dpc$LD1[grep( "AAAA", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "HHAA", row.names(dpc))]),
    y=1.1, yend=1.1,
    lty=1,
     size = 1, arrow = arrow(length = unit(0.2, "inches")),
    color="gray50")+
  geom_segment(
    x = mean(dpc$LD1[grep( "HHHH", row.names(dpc))]),
    xend = mean(dpc$LD1[grep( "AAHH", row.names(dpc))]),
    y=1.5, yend=1.5,
     size = 1, arrow = arrow(length = unit(0.2, "inches")),
    color="black")+
  theme(legend.position = "none")

pdf("~/Documents/UVM/Reciprocal_transplant/figures/dapc_SNP_gen.pdf", height = 6, width = 4)
ggarrange(a2,b2,c2, ncol=1, nrow=3, common.legend=TRUE)
dev.off()


#############
### stats:

dpc.all <- dpc1
dpc.all$Line <- substr(row.names(dpc.all),3,4)
dpc.all$Treatment <- substr(row.names(dpc.all),1,2)
dpc.all$gp <- c(rep(1.1,4), rep(1.5,8),rep(1.1,4))
dpc.all$LD1 <- dpc.all$LD1*-1

dpc.all$rep <- substr(row.names(dpc.all), 9,13)

dpc.all$home <- 1
dpc.all$home[grep("AAAA|HHHH", row.names(dpc.all))] <- 0

prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V=1, nu=0.002,alpha.mu=0, alpha.V=1000)))

mod1=MCMCglmm(LD1~Line+Line:home, random=~rep, prior= prior, data=dpc.all,nitt=75000, thin=25, burnin=5000)
summary(mod1)
#            post.mean l-95% CI u-95% CI eff.samp    pMCMC
#(Intercept)    4.2464   3.0894   5.4160     2800 0.000714 ***
#LineHH        -8.4911  -9.8337  -7.1083     2800  < 4e-04 ***
#LineAA:home   -2.2096  -3.6851  -0.9245     3014 0.005714 **
#LineHH:home    2.1664   0.8041   3.5429     2800 0.001429 **

posterior.mode(mod1$VCV)

HPDinterval(mod1$Sol)

#awayDelta=abs(mod1$Sol[,"LineAA:home"])-abs(mod1$Sol[,"LineHH:home"])
awayDelta=abs(mod1$Sol[,"LineHH:home"])-abs(mod1$Sol[,"LineAA:home"])

# 95% credible interval:
HPDinterval(awayDelta)

if (is.na(table(awayDelta<0)[2])) {
  cat("p <",signif(1/length(awayDelta),1))
    } else {
    cat("p =",signif(table(awayDelta<0)[2]/length(awayDelta),2))
  } 


# F2

dpc.all <- dpc2
dpc.all$Line <- substr(row.names(dpc.all),3,4)
dpc.all$Treatment <- substr(row.names(dpc.all),1,2)
dpc.all$gp <- c(rep(1.1,4), rep(1.5,8),rep(1.1,4))
dpc.all$LD1 <- dpc.all$LD1*-1

dpc.all$rep <- substr(row.names(dpc.all), 9,13)

dpc.all$home <- 1
dpc.all$home[grep("AAAA|HHHH", row.names(dpc.all))] <- 0

prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V=1, nu=0.002,alpha.mu=0, alpha.V=1000)))

mod2=MCMCglmm(LD1~Line+Line:home, random=~rep, prior= prior, data=dpc.all,nitt=75000, thin=25, burnin=5000)
summary(mod2)

#            post.mean l-95% CI u-95% CI eff.samp   pMCMC
# (Intercept)   -3.7183  -4.8065  -2.4644     2800 0.00143 **
# LineHH         7.4262   6.2932   8.5262     2513 < 4e-04 ***
# LineAA:home    1.8559   0.7622   2.9456     2800 0.00286 **
# LineHH:home   -2.1932  -3.3456  -1.0985     2800 0.00143 **



posterior.mode(mod2$VCV)

HPDinterval(mod2$Sol)

#awayDelta=abs(mod1$Sol[,"LineAA:home"])-abs(mod1$Sol[,"LineHH:home"])
awayDelta=abs(mod1$Sol[,"LineHH:home"])-abs(mod1$Sol[,"LineAA:home"])

# 95% credible interval:
HPDinterval(awayDelta)

if (is.na(table(awayDelta<0)[2])) {
  cat("p <",signif(1/length(awayDelta),1))
    } else {
    cat("p =",signif(table(awayDelta<0)[2]/length(awayDelta),2))
  } 




## F3

dpc.all <- dpc3
dpc.all$Line <- substr(row.names(dpc.all),3,4)
dpc.all$Treatment <- substr(row.names(dpc.all),1,2)
dpc.all$gp <- c(rep(1.1,4), rep(1.5,8),rep(1.1,4))
dpc.all$LD1 <- dpc.all$LD1*-1
dpc.all$rep <- substr(row.names(dpc.all), 9,13)

dpc.all$home <- 1
dpc.all$home[grep("AAAA|HHHH", row.names(dpc.all))] <- 0

prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V=1, nu=0.002,alpha.mu=0, alpha.V=1000)))

mod3=MCMCglmm(LD1~Line+Line:home, random=~rep, prior= prior, data=dpc.all,nitt=75000, thin=25, burnin=5000)
summary(mod3)
#            post.mean l-95% CI u-95% CI eff.samp    pMCMC
#(Intercept)  -2.99296 -4.12458 -1.80411     2800 0.002143 **
#LineHH        6.00740  4.90780  7.05361     2800  < 4e-04 ***
#LineAA:home   1.20025  0.07212  2.13389     3439 0.029286 *
#LineHH:home  -2.70971 -3.78055 -1.73574     2800 0.000714 ***
posterior.mode(mod3$VCV)

HPDinterval(mod3$Sol)

#awayDelta=abs(mod1$Sol[,"LineAA:home"])-abs(mod1$Sol[,"LineHH:home"])
awayDelta=abs(mod1$Sol[,"LineHH:home"])-abs(mod1$Sol[,"LineAA:home"])

# 95% credible interval:
HPDinterval(awayDelta)

if (is.na(table(awayDelta<0)[2])) {
  cat("p <",signif(1/length(awayDelta),1))
    } else {
    cat("p =",signif(table(awayDelta<0)[2]/length(awayDelta),2))
  }




# correct for multiple testing:
# f1 
#LineAA:home   -2.2096  -3.6851  -0.9245     3014 0.005714 **
#LineHH:home    2.1664   0.8041   3.5429     2800 0.001429 **

#f2
# LineAA:home    1.8559   0.7622   2.9456     2800 0.00286 **
# LineHH:home   -2.1932  -3.3456  -1.0985     2800 0.00143 **

#F3
#LineAA:home   1.20025  0.07212  2.13389     3439 0.029286 *
#LineHH:home  -2.70971 -3.78055 -1.73574     2800 0.000714 ***
p.adjust(c(0.005714,0.001429,0.00286,0.00143,0.029286,0.000714), method = "fdr")
# [1] 0.034284 0.008574 0.017160 0.008580 0.175716 0.004284
