dat <- read.table("~/reciprocal_t/analysis/snp_loadings.txt", header=TRUE)
dat$snp <- as.character(dat$snp)
Split <- (strsplit(dat$snp, "\\:"))

id <- sapply(Split, "[", 1)
snp <- sapply(Split, "[", 2)

#### mean var.contr per gene

gene <- unique(id)
out <- as.data.frame(matrix(nrow=length(gene), ncol=2))
colnames(out) <- c("gene", "mean_LD1", colnames(dat)[2:ncol(dat)])
for(i in 1:length(gene)){
    a <- dat[which(id == gene[i]),]
    out[i,1] <- gene[i]
    out[i,2] <- mean(a$LD1)

    if(i%%1000 ==0){print(i)}
}

write.table(file="~/reciprocal_t/analysis/GO_enrich/snp_meanLD_go.txt", out, col.names=TRUE,
    row.names=FALSE, quote=FALSE,sep=",")
