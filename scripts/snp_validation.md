# SNP validation

Compare snp calls from DNA and RNA to validate our RNA calls here. 

Using data from:

Dam, H.G. et al. 2021. Rapid, but limited, adaptation of marine zooplankton to simultaneous warming and acidification. Nature Climate Change. 11: 780–786



### first aligning DNA to transcriptome, to get everything on same footing.

```bash


nohup bash ~/reciprocal_t/scripts/align_transcriptome.sh 2> ~/reciprocal_t/log_out/align_transcriptome.stderr_$(date +"%F_%R").txt 1> ~/reciprocal_t/log_out/align_transcriptome.stdout_$(date +"%F_%R").txt &

echo $! > ~/reciprocal_t/log_out/align_transcriptome.pid


## merge bams

cd ~/reciprocal_t/analysis/SNPvalidation

for i in $(ls *.bam | cut -f 1 -d "." | uniq )

do {

  samtools merge ~/reciprocal_t/analysis/SNPvalidation/merged/${i}.bam\
      ~/reciprocal_t/analysis/SNPvalidation/${i}.bam \
      ~/reciprocal_t/analysis/SNPvalidation/${i}.lane2.bam

echo "done with {i}"

}
done


```



```r

### call snps:

#!/bin/bash -l

nohup bash ~/reciprocal_t/scripts/varscan_val.sh 2> ~/reciprocal_t/log_out/varscan_val.stderr_$(date +"%F_%R").txt 1> ~/reciprocal_t/log_out/varscan_val.stdout_$(date +"%F_%R").txt &

echo $! > ~/reciprocal_t/log_out/varscan_val.pid

1621135 variant positions reported (1621135 SNP, 10302 indel)

### filter snps
```


```r

## filtering raw variants

library(stringr)

dat <- read.table("~/reciprocal_t/analysis/SNPvalidation/snp_validation_out", stringsAsFactors=FALSE, skip=1)

datnm <- read.table("~/reciprocal_t/analysis/SNPvalidation/snp_validation_out", stringsAsFactors=FALSE, nrows=1)

pops <- c(
                 "AA_F25_REP1","AA_F25_REP2","AA_F25_REP3","AA_F25_REP4","HH_F25_REP1","HH_F25_REP2","HH_F25_REP3","HH_F25_REP4")

colnames(dat) <- c(datnm[1,1:10], pops)
nrow(dat)
#[1] 1621135
# first, remove all where the number of samples not covered/called is 0
dat1 <- dat[which(dat$SamplesNC == 0),]
nrow(dat1)
# [1] 487580

dat2 <- dat1

# filter by coverage:
filt <- (median(as.numeric((str_split_fixed(dat2[,5] , ":", n=6))[,2])))*3
filt <- 900*3
dat3 <- dat2[(which(as.numeric((str_split_fixed(dat2[,5] , ":", n=6))[,2]) < filt)),]
nrow(dat3)
#[1] 403578

#many of these are skewed by indels. only keep reads where depth of actual bialleleic snps > 40
# from the manual: Also, VarScan reports variants on a biallelic basis.
    #That is, for a given SNP call, the "reads1" column is the number of
    #reference-supporting reads (RD), and the "reads2" column is the number of
    #variant-supporting reads (AD).
    #There may be additional reads at that position showing other bases (SNP or indel variants).
    #If these other variants meet the calling criteria, they will be reported in
    #their own line. If not, it may look like you have "missing" reads.
# columns for each call are:
    #consensus genotype, total depth, num of read 1, num of read 2, allele freq, pvalue of snp.
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

low_cv <- (apply(keep, 1, function(x) {ifelse((length(which(x < 50)) > 0), FALSE, TRUE)}))
sum(low_cv)
#[1] 212575
dat4 <- dat3[low_cv,]
nrow(dat4)
dat3 <- dat4
nrow(dat3)
# [1] 212575

# to get rid of multiallele sites or deletions/insertions
dat4 <- dat3[(which(as.numeric(lapply(strsplit(dat3[,4],","), length)) == 1)),]
nrow(dat4)
# [1] 189067
dat3 <- dat4
# here calculate allele freqs
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

# get rid of invariant sites
dat4 <- dat3[(which(rowSums(af) > 0)),]
dat3 <- dat4
af <- af[(which(rowSums(af) > 0)),]
nrow(dat3)
#[1] 188720

af.out <- (cbind(paste(dat3$Chrom, dat3$Position, sep=":"),af))
af.out <- cbind(dat4$Chrom, dat4$Position,af.out)

colnames(af.out) <- c("CHR", "POS","SNP", colnames(af.out)[4:ncol(af.out)])

###########
######
# save filtered genotypes
#####
##########

write.table(file="~/reciprocal_t/analysis/SNPvalidation/filtered_allele_val_freqs.txt", af.out, sep="\t",
              row.names=FALSE, quote=FALSE)


```


### merge dna and rna snps

```r

library(scales)
library(poolSeq)
library(dplyr)

val <- read.table("~/reciprocal_t/analysis/SNPvalidation/filtered_allele_val_freqs.txt", sep="\t",
              header=T)


rna <- read.table("~/reciprocal_t/analysis/filtered_allele_freqs.txt", sep="\t",
              header=T)

dat <- merge(val, rna, by="SNP")
nrow(dat)
# 6938 snps

# figure out which ones are the CMH ones.

all_out <- read.table("~/reciprocal_t/analysis/AF_change.txt", header=TRUE)

dat <- merge(val, all_out, by="SNP")

dat$AA_F25_mean <- rowMeans(dat[,grep("AA_F25", colnames(dat))])
dat$HH_F25_mean <- rowMeans(dat[,grep("HH_F25", colnames(dat))])
sig.out <- dat[which(dat$sig_all == TRUE),]
nrow(sig.out)
# only 228 in this set


############################################
# plot difference between DNA and RNA estimates
############################################

#################
# simulate data

# 2 gens of drift,
# 50x coverage
# 20 indivs sampled
# pop size ~2000
# read in syn file:

# ne is 414, from tonsa genomics

#############
# simulate data
#############

# get starting allele freq:

startingAF <- dat$AA_F25_mean
effective_pop <- 414
samp_size <- 20
census_size <- 2000
F0_depth <- 75

nreps=2000

aa.out <- matrix(nrow=length(startingAF), ncol=nreps)

for(i in 1:nreps){
  # simulate the trajectory
  simTraj_1 <- wf.traj(p0=startingAF, Ne=effective_pop, t=c(2), s=0, h=0.5, haploid=FALSE, approximate=FALSE)

# then need to add noise due to sampling only a subset of individuals as well as subset available DNA (ie, coverage)
## first, due to sampling subset of indivs:

  simTraj_1a <- matrix(sample.alleles(simTraj_1, size=samp_size, mode="individuals", Ncensus=200), nrow=nrow(simTraj_1), dimnames=dimnames( simTraj_1))

  # introduce sampling variance to mimic Pool-seq of the entire population at coverage for each timepoint
  simTraj_1b <- cbind(simTraj_1,cbind((simTraj_1a),(sample.alleles(simTraj_1a[,1], size=75, mode="coverage"))))

  colnames(simTraj_1b) <- c("drift_AF", "Indiv_adj_AF", "Cov_adj_AF","Coverage")

    aa.out[,i] <- simTraj_1b$Cov_adj_AF

}



###### HH

startingAF <- dat$HH_F25_mean
effective_pop <- 414
samp_size <- 20
census_size <- 2000
F0_depth <- 75

nreps=2000

hh.out <- matrix(nrow=length(startingAF), ncol=nreps)

for(i in 1:nreps){
  # simulate the trajectory
  simTraj_1 <- wf.traj(p0=startingAF, Ne=effective_pop, t=c(2), s=0, h=0.5, haploid=FALSE, approximate=FALSE)

# then need to add noise due to sampling only a subset of individuals as well as subset available DNA (ie, coverage)
## first, due to sampling subset of indivs:

  simTraj_1a <- matrix(sample.alleles(simTraj_1, size=samp_size, mode="individuals", Ncensus=200), nrow=nrow(simTraj_1), dimnames=dimnames( simTraj_1))

  # introduce sampling variance to mimic Pool-seq of the entire population at coverage for each timepoint
  simTraj_1b <- cbind(simTraj_1,cbind((simTraj_1a),(sample.alleles(simTraj_1a[,1], size=75, mode="coverage"))))

  colnames(simTraj_1b) <- c("drift_AF", "Indiv_adj_AF", "Cov_adj_AF","Coverage")

    hh.out[,i] <- simTraj_1b$Cov_adj_AF

}

#####################
# for sig loci, do AF matched

x <- dat$AA_F25_mean
y <- dat$AAAA_F3_mean
sig.x <- sig.out$AA_F25_mean
sig.y <- sig.out$AAAA_F3_mean


hh.out.match <- matrix(nrow=nrow(sig.out), ncol=nreps)

# assign bin:
bin = cut(dat$AA_F25_mean, 
                    breaks=seq(from=0, to=1, by=0.02))

sig.bin <- bin[which(dat$sig_all == TRUE)]

for(i in 1:nreps){
  # simulate the trajectory

  tmp.x <- c()
  tmp.y <- c()

  for(j in 1:length(sig.bin)){
    tmp.idx <- sample(which(bin == sig.bin[j]), 1)
    tmp.x[j] <- x[tmp.idx]
    tmp.y[j] <- y[tmp.idx]
  }

    hh.out.match[,i] <- resid(lm(tmp.y~tmp.x))

    if(i%%100 == 0){print(i)}
}






##################################
########
# compare sim results to RNA results
# is variation similar?

# take diff from mean AA

aa.diff <- abs(aa.out - dat$AA_F25_mean)

# get CI


out <- list()
for (i in 1: ncol(aa.diff)){

    out[[i]] <- data.frame(x=density(aa.diff[,i], bw = 0.005)$x,
        y= density(aa.diff[,i], bw = 0.005)$y,
        bin = cut(density(aa.diff[,i], bw = 0.005)$x, 
                    breaks=seq(from=-0.3, to=0.6, by=0.005)))
}

out.new <- do.call(rbind, out)

densities.qtiles.aa <- out.new %>%
group_by(bin) %>%
  summarise(q05 = quantile(y, 0.025),
            q50 = quantile(y, 0.5),
            q95 = quantile(y, 0.975))

densities.qtiles.aa$x <- seq(from=-0.015, to=0.475, by=0.005)


#sim_mean_aa <- rowMeans(cbind(aa1$Cov_adj_AF, aa2$Cov_adj_AF,
#                       aa3$Cov_adj_AF, aa4$Cov_adj_AF))

# sig perm:

aa.diff.sig <- aa.diff[which(dat$sig_all == TRUE),]
out.sig <- list()
for (i in 1: ncol(aa.diff.sig)){

    out.sig[[i]] <- data.frame(x=density(aa.diff.sig[,i], bw = 0.005)$x,
        y= density(aa.diff.sig[,i], bw = 0.005)$y,
        bin = cut(density(aa.diff.sig[,i], bw = 0.005)$x, 
                    breaks=seq(from=-0.6, to=6, by=0.005)))
}

out.new <- do.call(rbind, out.sig)

densities.qtiles.sig.aa <- out.new %>%
group_by(bin) %>%
  summarise(q05 = quantile(y, 0.025),
            q50 = quantile(y, 0.5),
            q95 = quantile(y, 0.975))

densities.qtiles.sig.aa$x <- seq(from=-0.01, to=0.415, by=0.005)


#######
## hh

# take diff from mean

hh.diff <- abs(hh.out - dat$HH_F25_mean)

# get CI


out <- list()
for (i in 1: ncol(hh.diff)){

    out[[i]] <- data.frame(x=density(hh.diff[,i], bw = 0.005)$x,
        y= density(hh.diff[,i], bw = 0.005)$y,
        bin = cut(density(hh.diff[,i], bw = 0.005)$x, 
                    breaks=seq(from=-0.6, to=6, by=0.005)))
}

out.new <- do.call(rbind, out)

densities.qtiles.hh <- out.new %>%
group_by(bin) %>%
  summarise(q05 = quantile(y, 0.025),
            q50 = quantile(y, 0.5),
            q95 = quantile(y, 0.975))

densities.qtiles.hh$x <- seq(from=-0.01, to=0.45, by=0.005)


hh.diff.sig <- hh.diff[which(dat$sig_all == TRUE),]
out.sig <- list()
for (i in 1: ncol(hh.diff.sig)){

    out.sig[[i]] <- data.frame(x=density(hh.diff.sig[,i], bw = 0.005)$x,
        y= density(hh.diff.sig[,i], bw = 0.005)$y,
        bin = cut(density(hh.diff.sig[,i], bw = 0.005)$x, 
                    breaks=seq(from=-0.6, to=6, by=0.005)))
}

out.new <- do.call(rbind, out.sig)

densities.qtiles.sig.hh <- out.new %>%
group_by(bin) %>%
  summarise(q05 = quantile(y, 0.025),
            q50 = quantile(y, 0.5),
            q95 = quantile(y, 0.975))

densities.qtiles.sig.hh$x <- seq(from=-0.01, to=0.45, by=0.005)


# resid CIs


out.sig <- list()
for (i in 1: ncol(hh.out.match)){

    out.sig[[i]] <- data.frame(x=density(hh.out.match[,i], bw = 0.005)$x,
        y= density(hh.out.match[,i], bw = 0.005)$y,
        bin = cut(density(hh.out.match[,i], bw = 0.005)$x, 
                    breaks=seq(from=-1, to=1, by=0.005)))
}

out.new <- do.call(rbind, out.sig)

densities.qtiles.match.hh <- out.new %>%
group_by(bin) %>%
  summarise(q05 = quantile(y, 0.025),
            q50 = quantile(y, 0.5),
            q95 = quantile(y, 0.975))

densities.qtiles.match.hh$x <- seq(from=-0.905, to=0.72, by=0.005)






####################
### plot
####################


pdf("~/reciprocal_t/figures/dna_vs_rna_error_scatter.pdf", 
            height=10.5, width=7)

par(mfrow=c(3,2))

x <- dat$AA_F25_mean
y <- dat$AAAA_F3_mean
sig.x <- sig.out$AA_F25_mean
sig.y <- sig.out$AAAA_F3_mean

plot(x=x, y=y, col=alpha("black", 0.4), pch=19,
    ylab=c("RNA allele frequency"), xlab=c("DNA allele frequency"), 
    main=c(paste0("AM; rsq: ", round(summary(lm(y~x))$adj.r.squared, 3))))
points(x=sig.x, y=sig.y, col=alpha("firebrick3", 0.6), pch=19)
abline(lm(y~ x), col="dodgerblue3", lwd=2)
abline(lm(sig.y~ sig.x), col="firebrick3", lwd=2)



x <- dat$HH_F25_mean
y <- dat$HHHH_F3_mean
sig.x <- sig.out$HH_F25_mean
sig.y <- sig.out$HHHH_F3_mean

plot(x=x, y=y, col=alpha("black", 0.4), pch=19,
    ylab=c("RNA allele frequency"), xlab=c("DNA allele frequency"), 
    main=c(paste0("GH; rsq: ", round(summary(lm(y~x))$adj.r.squared, 3))))
points(x=sig.x, y=sig.y, col=alpha("firebrick3", 0.6), pch=19)
abline(lm(y~ x), col="dodgerblue3", lwd=2)
abline(lm(sig.y~ sig.x), col="firebrick3", lwd=2)


#################
# plot diffs in af 
#################

x <- dat$AA_F25_mean
y <- dat$AAAA_F3_mean
sig.x <- sig.out$AA_F25_mean
sig.y <- sig.out$AAAA_F3_mean

alldiff <- abs(x-y)
#simdiff <- x-aa2$Cov_adj_AF
sigdiff <- abs(sig.x-sig.y)
#simsig <- simdiff[which(dat$sig_all == TRUE)]

##### plot

# AM
hist(alldiff, breaks=30, col="grey60", main=c("AM"),
        xlab=c("DNA-RNA estimate"), xlim=c(-0.01,1), ylim=c(0, 21), freq=F)
#hist(simdiff, breaks=30, col=alpha("firebrick3", 0.5), 
#   xlab=c("DNA-RNA estimate"), xlim=c(-1,1), add=T, freq=F)

lines(density(alldiff), col="black",
        lwd=3)
#lines(density(sigdiff), col="firebrick3",
#       lwd=3)
polygon(x=c(densities.qtiles.aa$x,rev(densities.qtiles.aa$x)),
    y=c(densities.qtiles.aa$q05,rev(densities.qtiles.aa$q95)),
    col=alpha("dodgerblue3", alpha=0.7),border=NA)

lines(densities.qtiles.aa$x,densities.qtiles.aa$q50, col="dodgerblue3", lwd=2 )

legend(0.6, 15, legend=c("All loci", "Simulated"),
       col=c("black", "dodgerblue3"), lty=1, cex=0.8, lwd=3.5)

# HH


x <- dat$HH_F25_mean
y <- dat$HHHH_F3_mean
sig.x <- sig.out$HH_F25_mean
sig.y <- sig.out$HHHH_F3_mean

alldiff <- abs(x-y)
#simdiff <- x-aa2$Cov_adj_AF
sigdiff <- abs(sig.x-sig.y)

hist(alldiff, breaks=30, col="grey60", main=c("GH"),
        xlab=c("DNA-RNA estimate"), xlim=c(-0.01,1), ylim=c(0, 21), freq=F)
#hist(simdiff, breaks=30, col=alpha("firebrick3", 0.5), 
#   xlab=c("DNA-RNA estimate"), xlim=c(-1,1), add=T, freq=F)

lines(density(alldiff), col="black",
        lwd=3)

polygon(x=c(densities.qtiles.hh$x,rev(densities.qtiles.hh$x)),
    y=c(densities.qtiles.hh$q05,rev(densities.qtiles.hh$q95)),
    col=alpha("dodgerblue3", alpha=0.7),border=NA)

lines(densities.qtiles.hh$x,densities.qtiles.hh$q50, col="dodgerblue3", lwd=2 )




m1 <- lm(y~x)  #Create a linear model
m2 <- lm(sig.y ~ sig.x)
hist(resid(m1), breaks=30, col=NA, border=NA, freq=F, 
        ylim=c(0,12),
        xlim=c(-0.7, 0.7),
        xlab=c("Residual"),
        main=c("")) #A density plot

polygon(x=c(densities.qtiles.match.hh$x,rev(densities.qtiles.match.hh$x)),
    y=c(densities.qtiles.match.hh$q05,rev(densities.qtiles.match.hh$q95)),
    col=alpha("dodgerblue3", alpha=0.6),border=NA)

lines(density(resid(m1)), col="black",
        lwd=3)
lines(densities.qtiles.match.hh$x,densities.qtiles.match.hh$q50, 
    col="dodgerblue3", lwd=3 )
lines(density(resid(m2)), col="firebrick3",
        lwd=3)


dev.off()

```

