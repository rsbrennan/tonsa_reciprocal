
library(dplyr)
library(scales)
library(stringr)

# comparing DAPC contribs

dgeF1 <- read.csv("~/reciprocal_t/analysis/GO_enrich/dge_contrib_F1.go.in", header=TRUE, stringsAsFactors=FALSE)
dgeF2 <- read.csv("~/reciprocal_t/analysis/GO_enrich/dge_contrib_F2.go.in", header=TRUE, stringsAsFactors=FALSE)
dgeF3 <- read.csv("~/reciprocal_t/analysis/GO_enrich/dge_contrib_F3.go.in", header=TRUE, stringsAsFactors=FALSE)
colnames(dgeF1) <- c("gene", "LD1")
colnames(dgeF2) <- c("gene", "LD1")
colnames(dgeF3) <- c("gene", "LD1")

dgeF1F2 <- inner_join(dgeF1, dgeF2, by="gene")
dgeF2F3 <- inner_join(dgeF2, dgeF3, by="gene")
dgeF1F3 <- inner_join(dgeF1, dgeF3, by="gene")

png("~/reciprocal_t/figures/dge_gen_comp.png",width = 8, height = 3, res=300, units="in")

par(mfrow=c(1,3))

plot(dgeF1F2$LD1.x,dgeF1F2$LD1.y, pch=19, col=alpha("black", alpha=0.3),
    main="DGE: F1 vs. F2", xlab="F1 contributions", ylab="F2 contributions")
abline(lm(dgeF1F2$LD1.y ~dgeF1F2$LD1.x ), col="red")
text("cor=0.27", x= 0.005, y=0.003)

plot(dgeF2F3$LD1.x,dgeF2F3$LD1.y, pch=19, col=alpha("black", alpha=0.3),
    main="DGE: F2 vs. F3", xlab="F2 contributions", ylab="F3 contributions")
abline(lm(dgeF2F3$LD1.y ~dgeF2F3$LD1.x ), col="red")
text("cor=0.29", x= 0.005, y=0.003)

plot(dgeF1F3$LD1.x,dgeF1F3$LD1.y, pch=19, col=alpha("black", alpha=0.3),
main="DGE: F1 vs. F3", xlab="F1 contributions", ylab="F3 contributions")
abline(lm(dgeF1F3$LD1.y ~dgeF1F3$LD1.x ), col="red")
text("cor=0.40", x= 0.005, y=0.003)

dev.off()

cor.test(dgeF1F2$LD1.x, dgeF1F2$LD1.y)
cor.test(dgeF2F3$LD1.x, dgeF2F3$LD1.y)
cor.test(dgeF1F3$LD1.x, dgeF1F3$LD1.y)


##### snps

dat <- read.table("~/reciprocal_t/analysis/AF_change.txt", header=TRUE)

# some p-values are NA, bc no variation. change to 1. 
dat$pval_f1[which(is.na((dat$pval_f1)))] <- 1
dat$pval_f2[which(is.na((dat$pval_f2)))] <- 1
dat$pval_f3[which(is.na((dat$pval_f3)))] <- 1

dat$pval_f2[which(dat$pval_f2 == 0)] <- min(dat$pval_f2[which(dat$pval_f2 > 0)])
dat$pval_f3[which(dat$pval_f3 == 0)] <- min(dat$pval_f3[which(dat$pval_f3 > 0)])

dat$pval_f1[which(is.na((dat$pval_f1)))] <- 1
dat$pval_f2[which(is.na((dat$pval_f2)))] <- 1
dat$pval_f3[which(is.na((dat$pval_f3)))] <- 1

# look at overlap between snp and dge

F1 <- inner_join(dat, dgeF1, by="gene")
F2 <- inner_join(dat, dgeF2, by="gene")
F3 <- inner_join(dat, dgeF3, by="gene")

F1$pval_f1 <- -log10(F1$pval_f1)
F2$pval_f2 <- -log10(F2$pval_f2)
F3$pval_f3 <- -log10(F3$pval_f3)

#cor.test(F1$LD1, -log10(F1$pval_f1))

summary(lm(F1$pval_f1 ~ F1$LD1))
summary(lm(F2$pval_f2 ~ F2$LD1))
summary(lm(F3$pval_f3 ~ F3$LD1))

png("~/reciprocal_t/figures/snp_dge_comp.png",width = 8, height = 3, res=300, units="in")

par(mfrow=c(1,3))

plot(x= F1$LD1, y = (F1$pval_f1), pch=19, col=alpha("black", alpha=0.3),
    main="F1", xlab="DGE contributions", ylab="-log10 CMH p-value")
text("R-sq = 0.002", x= 0.002, y=200)

plot(F2$LD1, F2$pval_f2, pch=19, col=alpha("black", alpha=0.3),
    main="F2", xlab="DGE contributions", ylab=" ")
text("R-sq = 0.0007", x= 0.001, y=225)

plot(F3$LD1, F3$pval_f3, pch=19, col=alpha("black", alpha=0.3),
    main="F3", xlab="DGE contributions", ylab=" ")
text("R-sq = 0.003", x= 0.002, y=150)

dev.off()
