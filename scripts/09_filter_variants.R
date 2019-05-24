
## filtering raw variants

library(stringr)

dat <- read.table("~/reciprocal_t/analysis/snp_all_out", stringsAsFactors=FALSE, skip=1)

datnm <- read.table("~/reciprocal_t/analysis/snp_all_out", stringsAsFactors=FALSE, nrows=1)

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

colnames(dat) <- c(datnm[1,1:10], pops)
nrow(dat)
#[1] 5547802
# first, remove all where the number of samples not covered/called is 0
dat1 <- dat[which(dat$SamplesNC == 0),]
nrow(dat1)
# [1] 920486
# now, filter to include only variable sites
    # have 48 samples. want at least 4 variable.
    # so homozygotes both need to be < 44
#dat1 <- dat1[which(dat1$SamplesRef < 44 & dat1$SamplesHom < 44 ),]
#nrow(dat1)d
#[1] 682418d
dat2 <- dat1d

# filter by coverage:
filt <- (median(as.numeric((str_split_fixed(dat2[,5] , ":", n=6))[,2])))*3

dat3 <- dat2[(which(as.numeric((str_split_fixed(dat2[,5] , ":", n=6))[,2]) < filt)),]
nrow(dat3)
#[1] 797438

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
#[1] 451201
dat4 <- dat3[low_cv,]
nrow(dat4)
dat3 <- dat4
nrow(dat3)
# [1] 451201

# to get rid of multiallele sites or deletions/insertions
dat4 <- dat3[(which(as.numeric(lapply(strsplit(dat3[,4],","), length)) == 1)),]
nrow(dat4)
# [1] 411702
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
#[1] 411379

af.out <- (cbind(paste(dat3$Chrom, dat3$Position, sep=":"),af))

afct.maf <- (sapply(af,function(x)
          ifelse(x > 0.5, (1-x), x)))

low_maf <- (apply(afct.maf, 1, function(x) {ifelse((length(which(x > 0.025)) < 4), FALSE, TRUE)}))
sum(low_maf)
# [1] 322595

dat4 <- dat3[low_maf,]
nrow(dat4)
#[1] 322595

af_f <- af[low_maf,]

af.out <- (cbind(paste(dat4$Chrom, dat4$Position, sep=":"),af_f))

colnames(af.out) <- c("SNP", colnames(af_f))

###########
######
# save filtered genotypes
#####
##########

write.table(file="~/reciprocal_t/analysis/filtered_variants.txt", dat4, sep="\t",
              row.names=FALSE, quote=FALSE)

write.table(file="~/reciprocal_t/analysis/filtered_allele_freqs.txt", af.out, sep="\t",
              row.names=FALSE, quote=FALSE)

