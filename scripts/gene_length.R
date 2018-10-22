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
df <- read.table("/data/atonsa/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.gtf", header=FALSE)
nrow(df)

genes <- unique(dat3[,1])

length(genes)
#[1] 72511


out <- as.data.frame(matrix(ncol=2, nrow=length(genes)))
colnames(out) <- c("Chrom", "length")

for (i in 1:length(genes)){
    # subset to only gene of interest
    gene_sub <- df[which(df$V1 == genes[i]),]
    # add gene name and ref/var id, number of snps in gene
    out$Chrom[i] <- as.character(genes[i])
    # cycle through each population
    out$length[i] <- gene_sub$V5[nrow(gene_sub)] - gene_sub$V4[1]
    if(i%%1000 ==0){print(i)}
}

mean(out$length)
# [1] 1313.327

write.table(out,"~/reciprocal_t/analysis/gene_length.txt", sep="\t",
    quote=FALSE, row.names=FALSE)

hist(out$length, xlim=c(0,3000), breaks=100)

sum(out$length < 1000)

