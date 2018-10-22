
library(stringr)

dat <- read.table("~/reciprocal_t/analysis/snp_out", stringsAsFactors=FALSE, skip=1)

datnm <- read.table("~/reciprocal_t/analysis/snp_out", stringsAsFactors=FALSE, nrows=1)

colnames(dat) <- c(datnm[1,1:10], "AAAA_F1_REP1", "AAAA_F1_REP2", "AAAA_F1_REP3", "AAAA_F1_REP4","HHHH_F1_REP1", "HHHH_F1_REP2", "HHHH_F1_REP3", "HHHH_F1_REP4")

pops <- c("AAAA_F1_REP1", "AAAA_F1_REP2", "AAAA_F1_REP3", "AAAA_F1_REP4","HHHH_F1_REP1", "HHHH_F1_REP2", "HHHH_F1_REP3", "HHHH_F1_REP4")

# first, remove all where the number of samples not covered/called is not 8
dat1 <- dat[which(dat$SamplesNC == 0),]
nrow(dat1)
#[1] 708898

# now, filter to include only variable sites
dat1 <- dat1[which(dat1$SamplesRef < 6 & dat1$SamplesHom < 6 ),]
nrow(dat1)
#[1] 420010
dat2 <- dat1[which(dat1$SamplesHet > 0 | dat1$SamplesHom > 0),]
nrow(dat2)
# [1] 448070

# filter by coverage:
filt <- (median(as.numeric((str_split_fixed(dat2[,5] , ":", n=6))[,2])))*3

dat3 <- dat2[(which(as.numeric((str_split_fixed(dat2[,5] , ":", n=6))[,2]) < filt)),]
nrow(dat3)
#[1] 363200

#many of these are skewed by indels. only keep reads where depth of actual bialleleic snps > 40
# from the manual: Also, VarScan reports variants on a biallelic basis. That is, for a given SNP call, the "reads1" column is the number of reference-supporting reads (RD), and the "reads2" column is the number of variant-supporting reads (AD). There may be additional reads at that position showing other bases (SNP or indel variants). If these other variants meet the calling criteria, they will be reported in their own line. If not, it may look like you have "missing" reads. 
# columns for each call are: consensus genotype, total depth, num of read 1, num of read 2, allele freq, pvalue of snp.
keep <- as.data.frame(matrix(nrow=nrow(dat3), ncol=length(pops)))
colnames(keep) <- pops
    # cycle through each population
    for(i_pop in pops){

        tmp_pop <- dat3[,grep(i_pop, colnames(dat3))]
        maj <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,3])
        minor <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,4])

        # calculate af
        keep[,grep(i_pop, colnames(keep))] <- (maj+ minor)

    }

low_cv <- (apply(keep, 1, function(x) {ifelse((length(which(x < 30)) > 0), FALSE, TRUE)}))

sum(low_cv)

dat4 <- dat3[low_cv,]
nrow(dat4)
dat3 <- dat4
nrow(dat3)
# 362, 311

# columns for each call are: consensus genotype, total depth, num of read 1, num of read 2, allele freq, pvalue of snp.
af <- as.data.frame(matrix(nrow=nrow(dat3), ncol=length(pops)))
colnames(af) <- pops
    # cycle through each population
    for(i_pop in pops){

        tmp_pop <- dat3[,grep(i_pop, colnames(dat3))]
        maj <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,3])
        minor <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,4])

        # calculate af
        af[,grep(i_pop, colnames(af))] <- maj/(maj+ minor)

    }

#the first row should contain headers (a combi of locus and allele names)
#the first column should contain the names of the populations
freqs <- t(af)
colnames(freqs) <- c(paste(dat3$Chrom, dat3$Position, sep=":"))
nm <- cbind(colnames(af), freqs)
nm <- as.data.frame(cbind(substr(colnames(af), 1,4), nm)[,1:2])
colnames(nm) <- c("line", "id")

# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 
pcaResult <- prcomp(freqs, scale=TRUE) 
round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)
summary(pcaResult)

loadings <- as.data.frame(pcaResult$rotation)

hist(loadings$PC1, breaks=40)
qqnorm(loadings$PC1)
qqline(loadings$PC1, col = "red")

#This should show how the variables contribute to the PCA axes. e.g., negative values indicate negative relationship.

loads <- (cbind(row.names(loadings), loadings))
colnames(loads) <- c("gene", colnames(loadings)) 
write.table(loads, "~/reciprocal_t/analysis/pca_loadings.txt", sep="\t", quote=FALSE, row.names=FALSE)
loadings$PC1 <- abs(loadings$PC1)
# parse loadings to only include the most extreme value per gene
genes <- unique(str_split_fixed(row.names(loadings), ":", n=2)[,1])
loadings$gene <- str_split_fixed(row.names(loadings), ":", n=2)[,1]
new.loads <- as.data.frame(matrix(ncol=2, nrow=length(genes)))
colnames(new.loads) <- c("gene", "PC1")
for(i in 1:length(genes)){
    tmp_df <- loadings[which(loadings$gene == genes[i]),]
    new.loads$PC1[i] <- tmp_df$PC1[which.max(abs(tmp_df$PC1))]
    new.loads$gene[i] <- genes[i]
}

hist(new.loads$PC1)
write.table(new.loads, "~/reciprocal_t/analysis/pca_loadings_goenrich.txt", sep="\t", quote=FALSE, row.names=FALSE)

### use z scores #################

# with loadings$PC1
z.out <- as.data.frame(matrix(ncol=3, nrow=length(loadings$PC1)))
colnames(z.out) <- c("genes","z_score", "pval")
z.out$genes <- row.names(loadings)

tmp_mean <- mean(loadings$PC1)
tmp_sd <- sd(loadings$PC1)

z.out$z_score <- (sapply(loadings$PC1, function(x) (sum(x)- tmp_mean)/tmp_sd))
z.out$pval <- (sapply(z.out$z_score, function(x) 2*pnorm(-abs(x))))

hist(z.out$z_score, breaks=50, add=F, col="red", freq=F)
hist(loadings$PC1, breaks=40, freq=F)
hist(z.out$pval, breaks=50, add=F, col="red", freq=F)

# parse zscores to only include the most extreme value per gene
genes <- unique(str_split_fixed(z.out$genes, ":", n=2)[,1])
z.out$gene <- str_split_fixed(z.out$genes, ":", n=2)[,1]
new.z <- as.data.frame(matrix(ncol=3, nrow=length(genes)))
colnames(new.z) <- c("gene", "z_score", "pval")
for(i in 1:length(genes)){
    tmp_df <- z.out[which(z.out$gene == genes[i]),]
    new.z$z_score[i] <- tmp_df$z_score[which.max(abs(tmp_df$z_score))]
    new.z$pval[i] <- tmp_df$pval[which.min((tmp_df$pval))]
    new.z$gene[i] <- genes[i]
    if(i%%1000 == 0 ){print(i)}
}

par(mfrow=c(1,2))

hist(new.z$pval, breaks=50)
hist(new.z$z_score,breaks=45, col="red")
hist(new.loads$PC1,breaks=25, col="grey")

plot(new.z$z_score,new.loads$PC1 )
plot(z.out$z_score,loadings$PC1 )

all <- (merge(new.z, new.loads, by="gene"))
all[which(all$z_score < 0 & all$PC1 >0),]


write.table(new.z , "~/reciprocal_t/analysis/pca_zscores_goenrich.txt", sep="\t", quote=FALSE, row.names=FALSE)



#If you want the relative contribution of each variable, you can sum the total loadings for each PC axis (use absolute value for negatives) then divide each value with the column sum

prop_load <- as.data.frame(t(t(abs(loadings))/rowSums(t(abs(loadings))))*100)

# plot(cumsum(prop_load$PC1[order(prop_load$PC1, decreasing=TRUE)]))

####
##
## plot pca
##
####

library(ggplot2)
library(reshape)

# get proportion of total variance explained:
round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

dat.p <- data.frame(id=pops, Line=substr(pops, 1,4), PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,2])

a <- ggplot(dat.p, aes(PC1, PC2, fill=Line, shape=Line, size=Line)) +
        geom_point() +
        xlab(paste0("PC1: 23.9% variance")) +
        ylab(paste0("PC2: 17.4% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,24))+
        #scale_color_manual(values=c('brown3','blue3', 'darkgoldenrod3', "darkorchid3"))+
        scale_fill_manual(values=c('cornflowerblue', 'brown2'))+
        #theme(legend.position = c(0.88,0.17))+
        theme(legend.text=element_text(size=8),legend.title=element_blank())+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
        guides(fill=guide_legend(override.aes=list(shape=c(21,24), size=c(7,5), fill=c('cornflowerblue', 'brown2'))),
                shape=FALSE,
                size=FALSE)+
        scale_size_manual(values=c(7,5))
a

png("~/reciprocal_t/figures/pca.png", res=300, height=4, width=5, units="in")

a

dev.off()


# look at what loadings actually mean:

af[which.min(loads$PC1),]
af[which.max(loads$PC1),]

a1<- af[order(loads$PC1),]

a1<- t(a1[1:10,])
nm <- (as.character(substr(row.names(a1),1,4)))
nm <- cbind(nm, as.character(row.names(a1)))
adf <- as.data.frame(cbind(nm, a1))
colnames(adf) <- c("line", "id", colnames(a1))

mdata <- melt(adf, id=c("line", "id"))
mdata$value <- as.numeric(as.character(mdata$value))
mdata$variable <- as.factor(mdata$variable)


png("~/reciprocal_t/figures/loadings_af_neg.png",width = 6, height = 6, res=300, units="in")

ggplot(mdata, aes(x=line, y=value, color=line)) +
    theme_bw() +
    #stat_summary(fun.data="mean_se", mapping=aes(group=line), size=2.5,alpha = 0.5, show.legend=FALSE)+
    geom_jitter(width=0.1, size=3.5) +
    scale_color_manual(values=c("black", "firebrick3")) +
    ggtitle("negative loadings") +
    theme(axis.text.x=element_text(size=14),
        axis.title.x=element_blank()) +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~ variable, ncol=3)
dev.off()

png("~/reciprocal_t/figures/loadings_af_pos.png",width = 6, height = 6, res=300, units="in")

a1<- af[order(loads$PC1),]
a1<- t(a1[(nrow(a1)-10):nrow(a1),])
nm <- (as.character(substr(row.names(a1),1,4)))
nm <- cbind(nm, as.character(row.names(a1)))
adf <- as.data.frame(cbind(nm, a1))
colnames(adf) <- c("line", "id", colnames(a1))


mdata <- melt(adf, id=c("line", "id"))
mdata$value <- as.numeric(as.character(mdata$value))
mdata$variable <- as.factor(mdata$variable)
#png("~/reciprocal_t/figures/het.01_quantile.png",width = 6, height = 6, res=300, units="in")

ggplot(mdata, aes(x=line, y=value, color=line)) +
    theme_bw() +
    #stat_summary(fun.data="mean_se", mapping=aes(group=line), size=2.5,alpha = 0.5, show.legend=FALSE)+
    geom_jitter(width=0.1, size=3.5) +
    scale_color_manual(values=c("black", "firebrick3")) +
    ggtitle("positive loadings") +
    theme(axis.text.x=element_text(size=14),
        axis.title.x=element_blank()) +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~ variable, ncol=3)

dev.off()


########### het of extreme 


#############################################
#############################################
##
## calculate heterozygosity
##
#############################################
#############################################


af$id <- paste(dat3$Chrom, dat3$Position, sep=":")

out <- as.data.frame(matrix(ncol=8, nrow=nrow(af)))
colnames(out) <- colnames(af)[1:8]

for(i in 1:8){
    out[,i] <- sapply(af[,i], function(x) {(x)*(1-x)*2})
}

out <- cbind(af$id, out)
colnames(out) <- c("id", colnames(af)[1:8])

cutoff <- quantile(abs(loadings$PC1), c(0.999))
outlier <- row.names(loadings)[which(abs(loadings$PC1) > cutoff[1])]
out.plot <- which( as.character(out$id) %in% outlier)

sub_out <- out[out.plot,]

boxmelt <- melt(sub_out[2:ncol(sub_out)])

par(mar = (c(5, 5, 1, 0) + 0.1), oma = c(0, 0, 1, 1)) #
boxplot(boxmelt$value~boxmelt$variable, outpch = 19, outcex = 0.7,
#   levels=c("AAAA-1","AAAA-2", "AAAA-3", "AAAA-4", "HHHH-1", "HHHH-2", "HHHH-3", "HHHH-4")), 
#   ylim=c(70,150), ylab="", 
    yaxt="n", ylab="heterozygosity", frame.plot=F, 
    col=c("cornflowerblue","cornflowerblue","cornflowerblue", "cornflowerblue", 'brown2', 'brown2', 'brown2', 'brown2'), 
    axes=F,
    names=c("AAAA-1","AAAA-2", "AAAA-3", "AAAA-4", "HHHH-1", "HHHH-2", "HHHH-3", "HHHH-4"),
        cex.main=1.5, cex.axis=1.2,
        cex.lab=1.2,
        mgp=c(1.5,1,0))

axis(2 , cex.axis=1,mgp=c(0,0.6,0.1), las=1)
axis(1, cex.axis=1.2, labels=rep(" ", 8),
    at=seq(1,8,1))
axis(1, cex.axis=1.2, at=c(1,2,3), labels=rep(" ", 3))
text(x=seq(1,8,1), y=par()$usr[3]-0.08*(par()$usr[4]-par()$usr[3]),
labels=c("AAAA-1","AAAA-2", "AAAA-3", "AAAA-4", "HHHH-1", "HHHH-2", "HHHH-3", "HHHH-4"), srt=45, adj=1, xpd=TRUE, cex=1.2)

### do the same as above, but on a gene level


# need  sum of the reads of the major alleles for all SNPs in that window
# and  sum of the reads of the minor alleles for all SNPs in that window

# from https://www.nature.com/articles/s41559-018-0611-6#ref-CR113
    # Red fox genome assembly identifies genomic regions associated with tame and aggressive behaviours
calc_Hp <- function(numMajor, numMinor){
    Hp <- (2 * sum(numMajor) * sum(numMinor)) / ((sum(numMajor) + sum(numMinor))^2)

    return(Hp)
}

# columns for each call are: consensus genotype, total depth, num of read 1, num of read 2, allele freq, pvalue of snp.

genes <- unique(dat3$Chrom)
length(genes)
#[1] 11026

df <- dat3

out <- as.data.frame(matrix(ncol=12, nrow=length(genes)))
colnames(out) <- c(colnames(dat3[1]), colnames(dat3[3:4]), "num_variants", colnames(dat3)[11:18])

for (i in 1:length(genes)){
    # subset to only gene of interest
    gene_sub <- df[which(df$Chrom == genes[i]),]
    # add gene name and ref/var id, number of snps in gene
    out$Chrom[i] <- genes[i]
    out[i,2:3] <- gene_sub[1,3:4]
    out$num_variants[i] <- nrow(gene_sub)
    # cycle through each population
    for(i_pop in pops){

        tmp_pop <- gene_sub[,grep(i_pop, colnames(gene_sub))]
        maj <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,3])
        minor <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,4])

        # calculate heterozygosity
        out[i,grep(i_pop, colnames(out))] <- calc_Hp(maj, minor)

    }
    if(i%%1000==0){print(i)}
}

# isolate variable genes only. 
variable <- out[which(out$num_variants > 10),]
nrow(variable)
# [1] 7898


library(stringr)
cutoff <- quantile(abs(loadings$PC1), c(0.999))
outlier <- row.names(loadings)[which(abs(loadings$PC1) > cutoff[1])]

outlier_gene <- str_split_fixed(outlier, ":", n=2)[,1]

new <- out[which(variable$Chrom %in% outlier_gene),]

boxmelt <- melt(new[5:ncol(new)])

par(mar = (c(5, 5, 1, 0) + 0.1), oma = c(0, 0, 1, 1)) #
boxplot(boxmelt$value~boxmelt$variable, outpch = 19, outcex = 0.7,
#   levels=c("AAAA-1","AAAA-2", "AAAA-3", "AAAA-4", "HHHH-1", "HHHH-2", "HHHH-3", "HHHH-4")), 
#   ylim=c(70,150), ylab="", 
    yaxt="n", ylab="heterozygosity", frame.plot=F, 
    col=c("cornflowerblue","cornflowerblue","cornflowerblue", "cornflowerblue", 'brown2', 'brown2', 'brown2', 'brown2'), 
    axes=F,
    names=c("AAAA-1","AAAA-2", "AAAA-3", "AAAA-4", "HHHH-1", "HHHH-2", "HHHH-3", "HHHH-4"),
        cex.main=1.5, cex.axis=1.2,
        cex.lab=1.2,
        mgp=c(1.5,1,0))

axis(2 , cex.axis=1,mgp=c(0,0.6,0.1), las=1)
axis(1, cex.axis=1.2, labels=rep(" ", 8),
    at=seq(1,8,1))
axis(1, cex.axis=1.2, at=c(1,2,3), labels=rep(" ", 3))
text(x=seq(1,8,1), y=par()$usr[3]-0.08*(par()$usr[4]-par()$usr[3]),
labels=c("AAAA-1","AAAA-2", "AAAA-3", "AAAA-4", "HHHH-1", "HHHH-2", "HHHH-3", "HHHH-4"), srt=45, adj=1, xpd=TRUE, cex=1.2)

