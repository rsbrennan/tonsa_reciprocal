
## pca, etc for full data set, f1, f2, f3.

library(stringr) # yes
library(ggplot2)
library(reshape)
library(data.table)
library(gridExtra)
library(scales)
library(ggpubr)

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

##########################################################################
##########################################################################
# CMH
##########################################################################
##########################################################################

# need count data for this.
# get from dat3

A1 <- dat3[11:ncol(dat3)]
A1[] <- lapply(A1, function(x){as.numeric(str_split_fixed(x, ":", n=6)[,3]) })
A2 <- dat3[11:ncol(dat3)]
A2[] <- lapply(A2, function(x){as.numeric(str_split_fixed(x, ":", n=6)[,4]) })

#hist(as.numeric(str_split_fixed(dat3[,5], ":", n=6)[,2]))

#hist(as.numeric(str_split_fixed(dat3[,5], ":", n=6)[,2])[which(pvals_f1 < 0.05/(length(pvals_f1)*3))])


#which(pvals_f1 < 0.05/(length(pvals_f1)*3))

pvals_f1 <- c()
pvals_f2 <- c()
pvals_f3 <- c()

for(i in 1:nrow(A1)){

    #pull out snp
    sub_A1 <- A1[i,]
    sub_A2 <- A2[i,]
    #transform to data frame
    sub_A1 <- stack(sub_A1)
    sub_A2 <- stack(sub_A2)
    # add ID
    sub_A1$allele <- rep("ac1", nrow(sub_A1))
    sub_A2$allele <- rep("ac2", nrow(sub_A2))
    # add col ID
    colnames(sub_A1) <- c("count", "ind", "allele")
    colnames(sub_A2) <- c("count", "ind", "allele")
    # combine all
    sub_all <- rbind(sub_A1, sub_A2)

    # add ids
    sub_all$group <- substr(sub_all$ind, 1,4)
    sub_all$replicate <- substr(sub_all$ind, 9,12)
    sub_all$generation <- substr(sub_all$ind, 6,7)

    #only using control lines. remove other lines:
    sub_all <- sub_all[grep("AAHH|HHAA", sub_all$group, invert=TRUE),]

    # pull out each generation
    sub_f1 <- sub_all[grep("F1", sub_all$generation),]
    sub_f2 <- sub_all[grep("F2", sub_all$generation),]
    sub_f3 <- sub_all[grep("F3", sub_all$generation),]

    # table for cmh
    f1.xtabs = xtabs(count ~ allele + group + replicate,
                data=sub_f1)
    f2.xtabs = xtabs(count ~ allele + group + replicate,
                data=sub_f2)
    f3.xtabs = xtabs(count ~ allele + group + replicate,
                data=sub_f3)

    pvals_f1[i] <- mantelhaen.test(f1.xtabs)$p.value
    pvals_f2[i] <- mantelhaen.test(f2.xtabs)$p.value
    pvals_f3[i] <- mantelhaen.test(f3.xtabs)$p.value

    if (i%%10000 == 0){print(i)}

}

# pull out allele freqs from full data set
length(which(pvals_f1 < 0.05/(length(pvals_f1)*3) & pvals_f2 < 0.05/(length(pvals_f2)*3) & pvals_f3 < 0.05/(length(pvals_f3)*3)))
# 17720

length(which(pvals_f1 < 0.05/(length(pvals_f1)*3)))
#[1] 59823
length(which(pvals_f2 < 0.05/(length(pvals_f2)*3)))
#[1] 74140
length(which(pvals_f3 < 0.05/(length(pvals_f3)*3)))
#[1] 69864

# calc mean af for all loci

pop.id <- unique(substr(colnames(af)[2:length(colnames(af))], 1,7))

mean.out <- as.data.frame(matrix(nrow=nrow(af),ncol=length(pop.id)))

colnames(mean.out) <- pop.id

for(i in 1:length(pop.id)){

    tmp_pop <- af[,grep(pop.id[i], colnames(af))]
    mean.out[,i] <- apply(tmp_pop, 1, mean)
    print(i)

}

# combine pvalues for each gen, overlapping, af for each, mean af.
colnames(mean.out) <- paste(colnames(mean.out), "mean", sep="_")

all_out <- cbind(af, mean.out)

all_out$pval_f1 <- pvals_f1
all_out$pval_f2 <- pvals_f2
all_out$pval_f3 <- pvals_f3

# add sig val
all_out$sig_f1 <- FALSE
all_out$sig_f2 <- FALSE
all_out$sig_f3 <- FALSE
all_out$sig_all <- FALSE

all_out$sig_f1[which(pvals_f1 < 0.05/(length(pvals_f1)*3))] <- TRUE
all_out$sig_f2[which(pvals_f2 < 0.05/(length(pvals_f2)*3))] <- TRUE
all_out$sig_f3[which(pvals_f3 < 0.05/(length(pvals_f3)*3))] <- TRUE

all_out$sig_all[which(pvals_f1 < 0.05/(length(pvals_f1)*3) & pvals_f2 < 0.05/(length(pvals_f2)*3) & pvals_f3 < 0.05/(length(pvals_f3)*3))] <- TRUE


sum(all_out$sig_f1)
#[1] 59823
sum(all_out$sig_f2)
#[1] 74140
sum(all_out$sig_f3)
#[1] 69864
sum(all_out$sig_all)
#[1] 17720

write.table(all_out, 
    file="~/reciprocal_t/analysis/AF_change.txt", sep="\t", quote=FALSE, row.names=FALSE) 

