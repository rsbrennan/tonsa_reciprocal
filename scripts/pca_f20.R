
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

dat_t <- t(af)

library(pcadapt)

filename <- read.pcadapt(dat_t,type="pool")
x <- pcadapt(filename,K=7) # calc actual pca

# get proportion of total variance explained:
x$singular.values[1]/sum(x$singular.values)*100
x$singular.values[2]/sum(x$singular.values)*100

dat.p <- data.frame(id=pops, line=substr(pops, 1,4), PC1 = x$scores[,1],  PC2= x$scores[,2])

sp <- c(21, 22)
bg.col <-  c("gray48","firebrick3")

png("~/reciprocal_t/figures/pca.png", res=300, height=6, width=6, units="in")

par(mfrow = c(1, 1), mar=c(3, 3, 1.7, 1), mgp=c(3, 1, 0), las=0)
plot(y=dat.p$PC2, x=dat.p$PC1,
    pch=sp[as.factor(dat.p$line)],
    cex=1.6, col="black",
    bg=bg.col[as.factor(dat.p$line)],
    ylab="",
    xlab="",
    cex.lab=1.1, cex.axis=1,
    xaxt="n",yaxt="n"
)

axis(1, mgp=c(2, .5, 0), cex.axis=0.7) # second is tick mark labels
axis(2, mgp=c(2, .5, 0), cex.axis=0.7)

title(xlab="PC1: 18.7%", line=2, cex.lab=0.9)
title(ylab="PC2: 16.0%", line=2, cex.lab=0.9)

legend("bottomright", pch=c(21,22),
    pt.bg=c("gray31", "firebrick3") ,
    legend=c("AAAA", "HHHH" ),
    pt.cex=1.6, cex=0.8)

dev.off()
