
dat <- read.table("~/reciprocal_t/analysis/snp_contrib_F1.txt", header=TRUE)
dat$snp <- as.character(dat$gene)
Split <- (strsplit(dat$snp, "\\:"))

id <- sapply(Split, "[", 1)
snp <- sapply(Split, "[", 2)

gene <- unique(id)
out <- as.data.frame(matrix(nrow=length(gene), ncol=2))
colnames(out) <- c("gene", "mean_LD1")
for(i in 1:length(gene)){
    a <- dat[which(id == gene[i]),]
    out[i,1] <- gene[i]
    out[i,2] <- max(a$LD1)
    if(i%%1000 ==0){print(i)}
}

write.table(file="~/reciprocal_t/analysis/GO_enrich/snp_maxContrib_F1.txt", out, col.names=TRUE,
    row.names=FALSE, quote=FALSE,sep=",")

dat <- read.table("~/reciprocal_t/analysis/snp_contrib_F2.txt", header=TRUE)
dat$snp <- as.character(dat$gene)
Split <- (strsplit(dat$snp, "\\:"))

id <- sapply(Split, "[", 1)
snp <- sapply(Split, "[", 2)

gene <- unique(id)
out <- as.data.frame(matrix(nrow=length(gene), ncol=2))
colnames(out) <- c("gene", "mean_LD1")
for(i in 1:length(gene)){
    a <- dat[which(id == gene[i]),]
    out[i,1] <- gene[i]
    out[i,2] <- max(a$LD1)
    if(i%%1000 ==0){print(i)}
}

write.table(file="~/reciprocal_t/analysis/GO_enrich/snp_maxContrib_F2.txt", out, col.names=TRUE,
    row.names=FALSE, quote=FALSE,sep=",")


dat <- read.table("~/reciprocal_t/analysis/snp_contrib_F3.txt", header=TRUE)
dat$snp <- as.character(dat$SNP)
Split <- (strsplit(dat$snp, "\\:"))

id <- sapply(Split, "[", 1)
snp <- sapply(Split, "[", 2)

gene <- unique(id)
out <- as.data.frame(matrix(nrow=length(gene), ncol=2))
colnames(out) <- c("gene", "mean_LD1")
for(i in 1:length(gene)){
    a <- dat[which(id == gene[i]),]
    out[i,1] <- gene[i]
    out[i,2] <- max(a$LD1)
    if(i%%1000 ==0){print(i)}
}

write.table(file="~/reciprocal_t/analysis/GO_enrich/snp_maxContrib_F3.txt", out, col.names=TRUE,
    row.names=FALSE, quote=FALSE,sep=",")


#######
# from cmh

# find max for every gene

dat <- read.table("~/reciprocal_t/analysis/AF_change.txt", header=TRUE)
gene <- unique(dat$gene)

dat$pval_f1[which(is.na(dat$pval_f1))] <- 1
dat$pval_f2[which(is.na(dat$pval_f2))] <- 1
dat$pval_f3[which(is.na(dat$pval_f3))] <- 1

out <- as.data.frame(matrix(nrow=length(gene), ncol=4))
colnames(out) <- c("gene", "f1", "f2", "f3")
for(i in 1:length(gene)){
    a <- dat[which(dat$gene == gene[i]),]
    out$gene[i] <- as.character(gene[i])
    out$f1[i] <- min(a$pval_f1, na.rm=TRUE)
    out$f2[i] <- min(a$pval_f2, na.rm=TRUE)
    out$f3[i] <- min(a$pval_f3, na.rm=TRUE)
    if(i%%1000 ==0){print(i)}
}

# make rank based p-vals

f1 <- (out[order(out$f1),])
f1$rank <- seq(1,nrow(f1), by=1)
f1$rank_pval <- f1$rank/nrow(f1)

f2 <- (out[order(out$f2),])
f2$rank <- seq(1,nrow(f2), by=1)
f2$rank_pval <- f2$rank/nrow(f2)

f3 <- (out[order(out$f3),])
f3$rank <- seq(1,nrow(f3), by=1)
f3$rank_pval <- f3$rank/nrow(f3)

out1 <- rbind(f1[,c(1,6)], f2[,c(1,6)])
out2 <- rbind(out1, f3[,c(1,6)])
out2$rank_pval <- -log10(out2$rank_pval)

write.table(file="~/reciprocal_t/analysis/GO_enrich/snp_maxContrib_all.txt", 
          cbind(c(out$gene,out$gene,out$gene), -log10(c(out$f1, out$f2, out$f3))), 
          col.names=TRUE,
    row.names=FALSE, quote=FALSE,sep=",")

write.table(file="~/reciprocal_t/analysis/GO_enrich/snp_rankPval_all.txt", 
          out2,
          col.names=TRUE,
    row.names=FALSE, quote=FALSE,sep=",")


# save specific generations


f1 <- out[,c(1,2)]
f2 <- out[,c(1,3)]
f3 <- out[,c(1,4)]

f1$f1 <- -log10(f1$f1)
f2$f2 <- -log10(f2$f2)
f3$f3 <- -log10(f3$f3)


write.table(file="~/reciprocal_t/analysis/GO_enrich/snp_cmhF1_all.txt", 
          f1,
          col.names=TRUE,
    row.names=FALSE, quote=FALSE,sep=",")

write.table(file="~/reciprocal_t/analysis/GO_enrich/snp_cmhF2_all.txt", 
          f2,
          col.names=TRUE,
    row.names=FALSE, quote=FALSE,sep=",")

write.table(file="~/reciprocal_t/analysis/GO_enrich/snp_cmhF3_all.txt", 
          f3,
          col.names=TRUE,
    row.names=FALSE, quote=FALSE,sep=",")

hist(-log10(out$f1))
hist(-log10(out$f2), add=T, col="blue")
hist(-log10(out$f3),add=T, col="red")

plot(-log10(out$f1),-log10(out$f2))
plot(-log10(out$f2),-log10(out$f3))

