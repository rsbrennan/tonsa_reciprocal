############################################################
############################################################
###
### compare loadings
###
############################################################
############################################################
library(scales)

dat <- read.table("~/reciprocal_t/analysis/snp_loadings.txt", header=TRUE)
snp <- dat[,1:2]
snp$gene <- (str_split_fixed(snp$snp, ":", n=2)[,1])

dge <- read.table("~/reciprocal_t/analysis/dge_loadings.txt", header=TRUE)

df <- merge(dge, snp, by="gene")
colnames(df) <- c("gene", "dge_loading", "gene2", "snp_loading")

png("~/reciprocal_t/figures/loading_corr.png", height = 5, width = 5, units="in", res=300)
par(mfrow = c(1, 1), mar=c(4, 4, 1.7, 1), mgp=c(3, 1, 0), las=0)

plot((df$dge_loading), (df$snp_loading),
    xlab="DGE: DAPC transcript contributions",
    ylab="SNP: DAPC SNP contributions",
    pch=19, col=alpha("black", 0.2))
abline(lm(df$snp_loading~df$dge_loading), col="red", lwd=3)

dev.off()

cor.test(df$snp_loading,df$dge_loading, method="pearson")

summary(lm(df$snp_loading~df$dge_loading))


# look at quantiles

dge.99 <- df$gene[which(df$dge_loading > quantile(df$dge_loading, 0.99))]
af.99 <- df$gene[which(df$snp_loading > quantile(df$snp_loading, 0.99))]

length(intersect(dge.99, af.99))
length(unique(dge.99))

dge.9 <- df$gene[which(df$dge_loading > quantile(df$dge_loading, 0.9))]
af.9 <- df$gene[which(df$af_loading > quantile(df$af_loading, 0.9))]

length(unique(intersect(dge.9, af.9)))
