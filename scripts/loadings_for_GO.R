
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
    out[i,2] <- mean(a$LD1)
    if(i%%1000 ==0){print(i)}
}

write.table(file="~/reciprocal_t/analysis/GO_enrich/snp_meanContrib_F1.txt", out, col.names=TRUE,
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
    out[i,2] <- mean(a$LD1)
    if(i%%1000 ==0){print(i)}
}

write.table(file="~/reciprocal_t/analysis/GO_enrich/snp_meanContrib_F2.txt", out, col.names=TRUE,
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
    out[i,2] <- mean(a$LD1)
    if(i%%1000 ==0){print(i)}
}

write.table(file="~/reciprocal_t/analysis/GO_enrich/snp_meanContrib_F3.txt", out, col.names=TRUE,
    row.names=FALSE, quote=FALSE,sep=",")

