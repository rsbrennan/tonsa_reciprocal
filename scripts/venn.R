df <- read.table("~/reciprocal_t/analysis/AF_change.txt", header=TRUE)

sum(df$sig_f1)
sum(df$sig_f2)
sum(df$sig_f3)

f1 <- df$SNP[which(df$sig_f1 == TRUE & df$sig_f2 == TRUE)
length(which(df$sig_f1 == TRUE & df$sig_f3 == TRUE))
length(which(df$sig_f2 == TRUE & df$sig_f3 == TRUE))
length(which(df$sig_f1 == TRUE & df$sig_f2 == TRUE& df$sig_f3 == TRUE))

# make venn

f1 <- as.character(df$SNP[which(df$sig_f1 == TRUE)])
f2 <- as.character(df$SNP[which(df$sig_f2 == TRUE)])
f3 <- as.character(df$SNP[which(df$sig_f3 == TRUE)])


library(VennDiagram)

v1 <- venn.diagram(x= list(F1=f1, F2=f2, F3=f3),
    filename=NULL,fill=c("dodgerblue3", "firebrick3", "springgreen3"),
        cat.cex = 1.5, cex=1.2)

pdf("~/reciprocal_t/figures/venn.pdf", h=5, w=5)
grid.newpage()
grid.draw(v1)
dev.off()

