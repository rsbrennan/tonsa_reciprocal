library(stringr)

dat <- read.table("~/reciprocal_t/analysis/snp_out", stringsAsFactors=FALSE, skip=1)

datnm <- read.table("~/reciprocal_t/analysis/snp_out", stringsAsFactors=FALSE, nrows=1)

colnames(dat) <- c(datnm[1,1:10], "AAAA_F1_REP1", "AAAA_F1_REP2", "AAAA_F1_REP3", "AAAA_F1_REP4","HHHH_F1_REP1", "HHHH_F1_REP2", "HHHH_F1_REP3", "HHHH_F1_REP4")

pops <- c("AAAA_F1_REP1", "AAAA_F1_REP2", "AAAA_F1_REP3", "AAAA_F1_REP4","HHHH_F1_REP1", "HHHH_F1_REP2", "HHHH_F1_REP3", "HHHH_F1_REP4")

# can filter easily by the summary columns

    # SamplesRef    Number of samples called reference (wildtype)
    # SamplesHet    Number of samples called heterozygous-variant
    # SamplesHom    Number of samples called homozygous-variant
    # SamplesNC Number of samples not covered / not called

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

median(as.numeric((str_split_fixed(dat3[,5] , ":", n=6))[,2]))/8

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

        # sum up reads
        keep[,grep(i_pop, colnames(keep))] <- (maj+ minor)

    }

low_cv <- (apply(keep, 1, function(x) {ifelse((length(which(x < 30)) > 0), FALSE, TRUE)}))

sum(low_cv)

dat4 <- dat3[low_cv,]
nrow(dat4)
dat3 <- dat4
nrow(dat3)
# [1] 362311


#############################################
#############################################
##
## calculate heterozygosity
##
#############################################
#############################################

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

####### 
#
# save output
#
#######

write.table(variable, "~/reciprocal_t/analysis/heterozygosity_output.txt", sep="\t", quote=FALSE, row.name=FALSE)



# plot overall distribution
png("~/reciprocal_t/figures/heterozygosity.png", res=300, height=6, width=6, units="in")

par(mfrow = c(1, 1), mar=c(3, 3, 3, 1), mgp=c(3, 1, 0), las=0)

plot(density( variable$AAAA_F1_REP1), ylim=c(0,5.6), 
    lwd=2, xaxt="n",yaxt="n", main="")
lines(density(variable$AAAA_F1_REP2), lwd=2)
lines(density(variable$AAAA_F1_REP3), lwd=2)
lines(density(variable$AAAA_F1_REP4), lwd=2)
lines(density(variable$HHHH_F1_REP1), lwd=2, col="red")
lines(density(variable$HHHH_F1_REP2), lwd=2, col="red")
lines(density(variable$HHHH_F1_REP3), lwd=2, col="red")
lines(density(variable$HHHH_F1_REP4), lwd=2, col="red")
legend("topright", c("AAAA", "HHHH"), 
    col=c("black", "red"), lwd=2, cex=0.8 )

axis(1, mgp=c(2, .5, 0), cex.axis=0.7) # second is tick mark labels
axis(2, mgp=c(2, .5, 0), cex.axis=0.7)

title(xlab="heterozygosity", line=2, cex.lab=0.9)
title(ylab="density", line=2, cex.lab=0.9)
title(main="heterozygosity: window size = transcript", line=1, cex.lab=0.9)

dev.off()

mean(variable$AAAA_F1_REP1)
mean(variable$AAAA_F1_REP2)
mean(variable$AAAA_F1_REP3)
mean(variable$AAAA_F1_REP4)
mean(variable$HHHH_F1_REP1)
mean(variable$HHHH_F1_REP2)
mean(variable$HHHH_F1_REP3)
mean(variable$HHHH_F1_REP4)


# boxplots
library(reshape)

boxmelt <- melt(variable[5:ncol(variable)])

png("~/reciprocal_t/figures/het_boxplot.png", res=300, height=5, width=6, units="in")

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

dev.off()

png("~/reciprocal_t/figures/heterozygosity.png", res=300, height=6, width=6, units="in")

par(mfrow = c(1, 1), mar=c(3, 3, 3, 1), mgp=c(3, 1, 0), las=0)

library(ggplot2)
# Basic box plot
p <- ggplot(boxmelt, aes(x=variable, y=value, fill=line)) + 
  geom_boxplot() +
  geom_point(alpha=0.01, size = 5, shape = 21) +
  theme_bw()
p

boxplot(boxmelt$value ~ boxmelt$variable)

, ylim=c(0,5.6), 
    lwd=2, xaxt="n",yaxt="n", main="")
lines(density(variable$AAAA_F1_REP2), lwd=2)
lines(density(variable$AAAA_F1_REP3), lwd=2)
lines(density(variable$AAAA_F1_REP4), lwd=2)
lines(density(variable$HHHH_F1_REP1), lwd=2, col="red")
lines(density(variable$HHHH_F1_REP2), lwd=2, col="red")
lines(density(variable$HHHH_F1_REP3), lwd=2, col="red")
lines(density(variable$HHHH_F1_REP4), lwd=2, col="red")
legend("topright", c("AAAA", "HHHH"), 
    col=c("black", "red"), lwd=2, cex=0.8 )

axis(1, mgp=c(2, .5, 0), cex.axis=0.7) # second is tick mark labels
axis(2, mgp=c(2, .5, 0), cex.axis=0.7)

title(xlab="heterozygosity", line=2, cex.lab=0.9)
title(ylab="density", line=2, cex.lab=0.9)
title(main="heterozygosity: window size = transcript", line=1, cex.lab=0.9)

dev.off()


##### low het genes

a1 <- quantile(variable$AAAA_F1_REP1, na.rm=TRUE, c(0.01, 0.05))
a2 <- quantile(variable$AAAA_F1_REP2, na.rm=TRUE, c(0.01, 0.05))
a3 <- quantile(variable$AAAA_F1_REP3, na.rm=TRUE, c(0.01, 0.05))
a4 <- quantile(variable$AAAA_F1_REP4, na.rm=TRUE, c(0.01, 0.05))
h1 <- quantile(variable$HHHH_F1_REP1, na.rm=TRUE, c(0.01, 0.05))
h2 <- quantile(variable$HHHH_F1_REP2, na.rm=TRUE, c(0.01, 0.05))
h3 <- quantile(variable$HHHH_F1_REP3, na.rm=TRUE, c(0.01, 0.05))
h4 <- quantile(variable$HHHH_F1_REP4, na.rm=TRUE, c(0.01, 0.05))

gdat_05 <- data.frame(line=c(rep("AAAA", 4), rep("HHHH", 4)),
            pi=c(a1[2], a2[2], a3[2], a4[2], h1[2], h2[2], h3[2], h4[2]))

gdat <- data.frame(line=c(rep("AAAA", 4), rep("HHHH", 4)),
            pi=c(a1[1], a2[1], a3[1], a4[1], h1[1], h2[1], h3[1], h4[1]))

library(ggplot2)

png("~/reciprocal_t/figures/het.01_quantile.png",width = 6, height = 6, res=300, units="in")

ggplot(gdat, aes(x=line, y=pi, color=line)) +
    theme_bw() +
    stat_summary(fun.data="mean_se", mapping=aes(group=line), size=2.5,alpha = 0.5, show.legend=FALSE)+
    geom_jitter(width=0.1, size=3.5) +
    scale_color_manual(values=c("black", "firebrick3")) +
    ggtitle("heterozygosity: 0.01 quantile") +
    theme(axis.text.x=element_text(size=14),
        axis.title.x=element_blank()) +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5))

dev.off()

png("~/reciprocal_t/figures/het.05_quantile.png",width = 6, height = 6, res=300, units="in")
ggplot(gdat_05, aes(x=line, y=pi, color=line)) +
    theme_bw() +
    stat_summary(fun.data="mean_se", mapping=aes(group=line), size=2.5,alpha = 0.5, show.legend=FALSE)+
    geom_jitter(width=0.1, size=3.5) +
    scale_color_manual(values=c("black", "firebrick3")) +
    ggtitle("heterozygosity: 0.05 quantile") +
    theme(axis.text.x=element_text(size=14),
        axis.title.x=element_blank()) +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5))
dev.off()




HHHH_mean <- rowMeans(variable[,9:12])
AAAA_mean <- rowMeans(variable[,5:8])

cutoff <- mean(quantile(AAAA_mean, 0.05), quantile(HHHH_mean, 0.05))

AAAA_REP1 <- which(variable$AAAA_F1_REP1 < quantile(variable$AAAA_F1_REP1, na.rm=TRUE, c(0.05)))
AAAA_REP2 <- which(variable$AAAA_F1_REP2 < quantile(variable$AAAA_F1_REP2, na.rm=TRUE, c(0.05)))
AAAA_REP3 <- which(variable$AAAA_F1_REP3 < quantile(variable$AAAA_F1_REP3, na.rm=TRUE, c(0.05)))
AAAA_REP4 <- which(variable$AAAA_F1_REP4 < quantile(variable$AAAA_F1_REP4, na.rm=TRUE, c(0.05)))

HHHH_REP1 <- which(variable$HHHH_F1_REP1 < quantile(variable$HHHH_F1_REP1, na.rm=TRUE, c(0.05)))
HHHH_REP2 <- which(variable$HHHH_F1_REP2 < quantile(variable$HHHH_F1_REP1, na.rm=TRUE, c(0.05)))
HHHH_REP3 <- which(variable$HHHH_F1_REP3 < quantile(variable$HHHH_F1_REP1, na.rm=TRUE, c(0.05)))
HHHH_REP4 <- which(variable$HHHH_F1_REP4 < quantile(variable$HHHH_F1_REP1, na.rm=TRUE, c(0.05)))

all_overlap <- (intersect(Reduce(intersect, list(HHHH_REP1, HHHH_REP2, HHHH_REP3, HHHH_REP4)), 
    Reduce(intersect, list(AAAA_REP1, AAAA_REP2, AAAA_REP3, AAAA_REP4))))
length(all_overlap)
# 87


##
## look for number of low het genes overlapping between 2, 3, and 4 for each group
## 

##### 2 reps
two.1 <- intersect(AAAA_REP1, AAAA_REP2)
two.2 <- intersect(AAAA_REP1, AAAA_REP3)
two.3 <- intersect(AAAA_REP1, AAAA_REP4)
two.4 <- intersect(AAAA_REP2, AAAA_REP3)
two.5 <- intersect(AAAA_REP2, AAAA_REP4)
two.6 <- intersect(AAAA_REP3, AAAA_REP4)
length(unique(c(two.1, two.2, two.3, two.4, two.5, two.6)))
# 415

two.1 <- intersect(HHHH_REP1, HHHH_REP2)
two.2 <- intersect(HHHH_REP1, HHHH_REP3)
two.3 <- intersect(HHHH_REP1, HHHH_REP4)
two.4 <- intersect(HHHH_REP2, HHHH_REP3)
two.5 <- intersect(HHHH_REP2, HHHH_REP4)
two.6 <- intersect(HHHH_REP3, HHHH_REP4)
length(unique(c(two.1, two.2, two.3, two.4, two.5, two.6)))
# 432

##### 3 reps
combinations <- as.data.frame(combn(c("AAAA_REP1", "AAAA_REP2", "AAAA_REP3", "AAAA_REP4"),3, simply=FALSE))


three.1 <- Reduce(intersect, list(AAAA_REP1, AAAA_REP2, AAAA_REP3))
three.2 <- Reduce(intersect, list(AAAA_REP1, AAAA_REP2, AAAA_REP4))
three.3 <- Reduce(intersect, list(AAAA_REP1, AAAA_REP3, AAAA_REP4))
three.4 <- Reduce(intersect, list(AAAA_REP4, AAAA_REP2, AAAA_REP3))
length(unique(c(three.1, three.2, three.3, three.4)))
# 253

three.1 <- Reduce(intersect, list(HHHH_REP1, HHHH_REP2, HHHH_REP3))
three.2 <- Reduce(intersect, list(HHHH_REP1, HHHH_REP2, HHHH_REP4))
three.3 <- Reduce(intersect, list(HHHH_REP1, HHHH_REP3, HHHH_REP4))
three.4 <- Reduce(intersect, list(HHHH_REP4, HHHH_REP2, HHHH_REP3))
length(unique(c(three.1, three.2, three.3, three.4)))
# 271

##### 4 reps

AAAA_overlap <- variable[(Reduce(intersect, list(AAAA_REP1, AAAA_REP2, AAAA_REP3, AAAA_REP4))),]
HHHH_overlap <- variable[(Reduce(intersect, list(HHHH_REP1, HHHH_REP2, HHHH_REP3, HHHH_REP4))),]
nrow(AAAA_overlap)
# 138
nrow(HHHH_overlap)
# 146


### determine outliers by permutation:

    # SNP positions are taken as fixed and allele frequencies are randomly shuffled across positions in each iteration
    # computing pooled heterozygosities for creeping windows of 40 K from the shuffled data and the genome wide lowest HP from a window ≥30 K formed by ≥10 SNPs is stored. After repeating this procedure for n iterations, the empirical threshold pertaining to error probability P = 0.001 is the value cutting of the 0.001 quantile in the ordered vector of minima.
    # from http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0049525#s4
    # and https://www.nature.com/articles/s41559-018-0611-6#ref-CR114


# subset df to only genes that ha

nrep = 10000

perm_out <- as.data.frame(matrix(ncol=8, nrow=nrep))
colnames(perm_out) <- pops

# df contains the read counts

df_var <- df[df$Chrom %in% variable$Chrom,]
genes <- unique(df_var$Chrom)


for(n_rep in 1:nrep){
    perm_df <- df_var

    for (pop_id in pops){
        perm_df[,grep(pop_id, colnames(perm_df))] <-  sample(perm_df[,grep(pop_id, colnames(perm_df))])
    }
    
    out <- as.data.frame(matrix(ncol=12, nrow=length(genes)))
    colnames(out) <- c(colnames(dat3[1]), colnames(dat3[3:4]), "num_variants", colnames(dat3)[11:18])
    
    for (i in 1:length(genes)){
        # subset to only gene of interest
        gene_sub <- perm_df[which(perm_df$Chrom == genes[i]),]
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

    #find lowest and add to output
    for(min_pop in pops){

        # calculate heterozygosity
        perm_out[n_rep,grep(min_pop, colnames(perm_out))] <- min(out[,grep(min_pop, colnames(out))] )
    
        }   

    print(paste("rep", n_rep, "done", sep=" "))
}





overlap_rm <- variable[-all_overlap,]

quantile((overlap_rm$AAAA_F1_REP1), c(0.01, 0.1, 0.2))
quantile((overlap_rm$AAAA_F1_REP2), c(0.01, 0.1, 0.2))
quantile((overlap_rm$AAAA_F1_REP3), c(0.01, 0.1, 0.2))
quantile((overlap_rm$AAAA_F1_REP4), c(0.01, 0.1, 0.2))
quantile((overlap_rm$HHHH_F1_REP1), c(0.01, 0.1, 0.2))
quantile((overlap_rm$HHHH_F1_REP2), c(0.01, 0.1, 0.2))
quantile((overlap_rm$HHHH_F1_REP3), c(0.01, 0.1, 0.2))
quantile((overlap_rm$HHHH_F1_REP4), c(0.01, 0.1, 0.2))



# z transform:

# ZHp = (Hp - μHp)/σHp

z.out <- as.data.frame(matrix(ncol=12, nrow=length(genes)))
z.out[,1:4] <- out[,1:4]
colnames(z.out) <- c(colnames(dat3[1]), colnames(dat3[3:4]), "num_variants", colnames(dat3)[11:18])

for(i_pop in pops){

    tmp_pop <- out[,grep(i_pop, colnames(out))]
    tmp_mean <- mean(tmp_pop)
    tmp_sd <- sd(tmp_pop)
    for(i in 1:length(tmp_pop)){
        z.out[i,grep(i_pop, colnames(z.out))] <- (tmp_pop[i] - tmp_mean)/tmp_sd
    }
}

