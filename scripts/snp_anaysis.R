
## pca, etc for full data set, f1, f2, f3.

library(stringr) # yes
library(ggplot2)
library(reshape)
library(data.table)
library(gridExtra)
library(scales)

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

# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 

pcaResult <- prcomp(freqs, scale=TRUE) 
#round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)
#summary(pcaResult)

####
##
## plot pca
##
####

# get proportion of total variance explained:
percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

dat.p <- data.frame(id=pops, Line=substr(pops, 1,4), 
        gen=substr(pops, 6,7),
        PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,2])

a <- ggplot(dat.p, aes(PC1, PC2, fill=Line, shape=gen)) +
        geom_point(size=4, color="black") +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,22,24))+
        #scale_color_manual(values=c('brown3','blue3', 'darkgoldenrod3', "darkorchid3"))+
        scale_fill_manual(values=c('brown3','blue3', 'darkgoldenrod3', "darkorchid3"))+
        #theme(legend.position = c(0.88,0.17))+
       # theme(legend.text=element_text(size=8),legend.title=element_blank())+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
       guides(fill=guide_legend(override.aes=list(shape=c(21,21,21, 21), 
                size=c(5,5,5,5), fill=c('brown3','blue3', 'darkgoldenrod3', "darkorchid3"))))+
       ggtitle("All generations")
        #        shape=FALSE,
        #        size=FALSE)+
        #scale_size_manual(values=c(7,5))
a


png("~/reciprocal_t/figures/pca_all.png", res=300, height=4, width=5, units="in")

a

dev.off()

# pc1 and 3

dat.p <- data.frame(id=pops, Line=substr(pops, 1,4), 
        gen=substr(pops, 6,7),
        PC1 = pcaResult$x[,1],  PC3= pcaResult$x[,3])

a2 <- ggplot(dat.p, aes(PC1, PC3, fill=Line, shape=gen)) +
        geom_point(size=4, color="black") +
        xlab(paste0("PC1: 10.6% variance")) +
        ylab(paste0("PC3: 4.5% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,22,24))+
        #scale_color_manual(values=c('brown3','blue3', 'darkgoldenrod3', "darkorchid3"))+
        scale_fill_manual(values=c('brown3','blue3', 'darkgoldenrod3', "darkorchid3"))+
        #theme(legend.position = c(0.88,0.17))+
       # theme(legend.text=element_text(size=8),legend.title=element_blank())+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
       guides(fill=guide_legend(override.aes=list(shape=c(21,21,21, 21), 
                size=c(5,5,5,5), fill=c('brown3','blue3', 'darkgoldenrod3', "darkorchid3"))))
        #        shape=FALSE,
        #        size=FALSE)+
        #scale_size_manual(values=c(7,5))
a2

png("~/reciprocal_t/figures/pca_all_pc1_3.png", res=300, height=4, width=5, units="in")

a

dev.off()


# f1 alone

freqs_f1 <- freqs[grep("F1", row.names(freqs)),]
freqs_f1 <- freqs_f1[,which(apply(freqs_f1, 2, var)>0)]
pcaResult <- prcomp(freqs_f1, scale=TRUE) 

percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

f1 <- pops[grep("F1", pops)]

dat.p <- data.frame(id=f1, Line=substr(f1, 3,4), 
    Treatment=substr(f1, 1,2),
        gen=substr(f1, 6,7),
        PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,2])

b <- ggplot(dat.p, aes(PC1, PC2, fill=Treatment, shape=Line, size=Line)) +
        geom_point() +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,24))+
        #scale_color_manual(values=c('brown3','blue3', 'darkgoldenrod3', "darkorchid3"))+
        scale_fill_manual(values=c('cornflowerblue', 'brown2'))+
        #theme(legend.position = c(0.88,0.17))+
        theme(legend.text=element_text(size=8))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
        ggtitle("F1")+
        guides(fill=guide_legend(override.aes=list(shape=21, size=7, fill=c('cornflowerblue', 'brown2')),order = 2),
                shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)),order = 1),
                size=FALSE)+
        scale_size_manual(values=c(7,5))
#b

# f2 alone

freqs_f2 <- freqs[grep("F2", row.names(freqs)),]
freqs_f2 <- freqs_f2[,which(apply(freqs_f2, 2, var)>0)]
pcaResult <- prcomp(freqs_f2, scale=TRUE) 

percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

f2 <- pops[grep("F2", pops)]

dat.p <- data.frame(id=f2, Line=substr(f2, 3,4), 
    Treatment=substr(f2, 1,2),
        gen=substr(f2, 6,7),
        PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,2])

c <- ggplot(dat.p, aes(PC1, PC2, fill=Treatment, shape=Line, size=Line)) +
        geom_point() +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,24))+
        #scale_color_manual(values=c('brown3','blue3', 'darkgoldenrod3', "darkorchid3"))+
        scale_fill_manual(values=c('cornflowerblue', 'brown2'))+
        #theme(legend.position = c(0.88,0.17))+
        theme(legend.text=element_text(size=8))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
        ggtitle("F2")+
        guides(fill=guide_legend(override.aes=list(shape=21, size=7, fill=c('cornflowerblue', 'brown2')),order = 2),
                shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)),order = 1),
                size=FALSE)+
        scale_size_manual(values=c(7,5))

#c

# f3 alone

freqs_f3 <- freqs[grep("F3", row.names(freqs)),]
freqs_f3 <- freqs_f3[,which(apply(freqs_f3, 2, var)>0)]
pcaResult <- prcomp(freqs_f3, scale=TRUE) 

percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

f3 <- pops[grep("F3", pops)]

dat.p <- data.frame(id=f3, Line=substr(f3, 3,4), 
    Treatment=substr(f3, 1,2),
        gen=substr(f3, 6,7),
        PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,2])

d <- ggplot(dat.p, aes(PC1, PC2, fill=Treatment, shape=Line, size=Line)) +
        geom_point() +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,24))+
        #scale_color_manual(values=c('brown3','blue3', 'darkgoldenrod3', "darkorchid3"))+
        scale_fill_manual(values=c('cornflowerblue', 'brown2'))+
        #theme(legend.position = c(0.88,0.17))+
        theme(legend.text=element_text(size=8))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
        ggtitle("F3")+
        guides(fill=guide_legend(override.aes=list(shape=21, size=7, fill=c('cornflowerblue', 'brown2')),order = 2),
                shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)),order = 1),
                size=FALSE)+
        scale_size_manual(values=c(7,5))
#d

png("~/reciprocal_t/figures/pca_all.png", res=300, height=6, width=12, units="in")

grid.arrange(a, a2, b, c, d, nrow = 2)

dev.off()


png("~/reciprocal_t/figures/pca_1_3.png", res=300, height=4, width=11, units="in")

ggarrange(b, c, d, nrow = 1, ncol=3, common.legend=TRUE)

dev.off()


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

hist(as.numeric(str_split_fixed(dat3[,5], ":", n=6)[,2]))

hist(as.numeric(str_split_fixed(dat3[,5], ":", n=6)[,2])[which(pvals_f1 < 0.05/(length(pvals_f1)*3))])


which(pvals_f1 < 0.05/(length(pvals_f1)*3))

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
af.out_f1 <- af[which(pvals_f1 < 0.05/(length(pvals_f1)*3)),]
af.out_f2 <- af[which(pvals_f2 < 0.05/(length(pvals_f2)*3)),]
af.out_f3 <- af[which(pvals_f3 < 0.05/(length(pvals_f3)*3)),]

length(which(pvals_f1 < 0.05/(length(pvals_f1)*3)))
#[1] 73272
length(which(pvals_f2 < 0.05/(length(pvals_f2)*3)))
#[1] 92575
length(which(pvals_f3 < 0.05/(length(pvals_f3)*3)))
#[1] 86926


# and plot these allele freqs
pop.id <- unique(substr(colnames(af.out)[2:length(colnames(af.out))], 1,7))

mean.out_f1 <- as.data.frame(matrix(nrow=nrow(af.out_f1),ncol=length(pop.id)))
mean.out_f2 <- as.data.frame(matrix(nrow=nrow(af.out_f2),ncol=length(pop.id)))
mean.out_f3 <- as.data.frame(matrix(nrow=nrow(af.out_f3),ncol=length(pop.id)))
colnames(mean.out_f1) <- pop.id
colnames(mean.out_f2) <- pop.id
colnames(mean.out_f3) <- pop.id

for(i in 1:length(pop.id)){

    tmp_pop <- af.out_f1[,grep(pop.id[i], colnames(af.out_f1))]
    mean.out_f1[,i] <- apply(tmp_pop, 1, mean)

    tmp_pop <- af.out_f2[,grep(pop.id[i], colnames(af.out_f2))]
    mean.out_f2[,i] <- apply(tmp_pop, 1, mean)

    tmp_pop <- af.out_f3[,grep(pop.id[i], colnames(af.out_f3))]
    mean.out_f3[,i] <- apply(tmp_pop, 1, mean)
}

#####
# plot
####

png("~/reciprocal_t/figures/trt_allele_freq.png", height = 7, width = 9, res=300, units="in")

par(mfrow=c(2,3))
x <- mean.out_f1$AAAA_F1 - mean.out_f1$HHHH_F1
y <- mean.out_f1$AAHH_F1 - mean.out_f1$HHHH_F1
plot(x =x,
    y = y, 
    col=alpha("black", 0.1), bg=alpha("black", 0.05), pch=21,
    ylim=c(-1,1), xlim=c(-1, 1),
    xlab="AAAA - HHHH",
    ylab="(HH in AA) - HHHH",
    main="F1")
abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="red", lwd=2)
text(x=-0.7, y=0.7, paste("Corr = ", 
  round(cor.test(x,y, method="pearson")$estimate, 3)))

x <- mean.out_f2$AAAA_F2 - mean.out_f2$HHHH_F2
y <- mean.out_f2$AAHH_F2 - mean.out_f2$HHHH_F2
plot(x = x,
    y = y, 
    col=alpha("black", 0.1), bg=alpha("black", 0.05), pch=21,
    ylim=c(-1,1), xlim=c(-1, 1),
    xlab="AAAA - HHHH",
    ylab="(HH in AA) - HHHH",
    main="F2")
abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="red", lwd=2)
text(x=-0.7, y=0.7, paste("Corr = ", 
  round(cor.test(x,y, method="pearson")$estimate, 3)))

x <- mean.out_f3$AAAA_F3 - mean.out_f3$HHHH_F3
y <- mean.out_f3$AAHH_F3 - mean.out_f3$HHHH_F3
plot(x = x,
    y = y, 
    col=alpha("black", 0.1), bg=alpha("black", 0.05), pch=21,
    ylim=c(-1,1), xlim=c(-1, 1),
    xlab="AAAA - HHHH",
    ylab="(HH in AA) - HHHH",
    main="F3")
abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="red", lwd=2)
text(x=-0.7, y=0.7, paste("Corr = ", 
  round(cor.test(x,y, method="pearson")$estimate, 3)))

# and for AA in HH
x <- mean.out_f1$HHHH_F1 - mean.out_f1$AAAA_F1
y <- mean.out_f1$HHAA_F1 - mean.out_f1$AAAA_F1
plot(x =x,
    y = y, 
    col=alpha("black", 0.1), bg=alpha("black", 0.05), pch=21,
    ylim=c(-1,1), xlim=c(-1, 1),
    xlab="HHHH - AAAA",
    ylab="(AA in HH) - AAAA",
    main="F1")
abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="red", lwd=2)
text(x=-0.7, y=0.7, paste("Corr = ", 
  round(cor.test(x,y, method="pearson")$estimate, 3)))

x <- mean.out_f2$HHHH_F2 - mean.out_f2$AAAA_F2
y <- mean.out_f2$HHAA_F2 - mean.out_f2$AAAA_F2
plot(x = x,
    y = y, 
    col=alpha("black", 0.1), bg=alpha("black", 0.05), pch=21,
    ylim=c(-1,1), xlim=c(-1, 1),
    xlab="HHHH - AAAA",
    ylab="(AA in HH) - AAAA",
    main="F2")
abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="red", lwd=2)
text(x=-0.7, y=0.7, paste("Corr = ", 
  round(cor.test(x,y, method="pearson")$estimate, 3)))

x <- mean.out_f3$HHHH_F3 - mean.out_f3$AAAA_F3
y <- mean.out_f3$HHAA_F3 - mean.out_f3$AAAA_F3
plot(x = x,
    y = y, 
    col=alpha("black", 0.1), bg=alpha("black", 0.05), pch=21,
    ylim=c(-1,1), xlim=c(-1, 1),
    xlab="HHHH - AAAA",
    ylab="(AA in HH) - AAAA",
    main="F3")
abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="red", lwd=2)

text(x=-0.7, y=0.7, paste("Corr = ", 
  round(cor.test(x,y, method="pearson")$estimate, 3)))

dev.off()




##########################################################################
##########################################################################
# loadings
##########################################################################
##########################################################################

freqs <- t(af)
colnames(freqs) <- c(paste(dat3$Chrom, dat3$Position, sep=":"))
nm <- cbind(colnames(af), freqs)
nm <- as.data.frame(cbind(substr(colnames(af), 1,4), nm)[,1:2])
colnames(nm) <- c("line", "id")

# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 

freqs.ctr <- freqs[grep("AAAA|HHHH", row.names(freqs)),]

freqs.ctr <- freqs.ctr[,which(apply(freqs.ctr, 2, var)>0)]
freqs.ctr[] <- lapply(freqs.ctr, function(x){as.numeric(freqs.ctr) })

pcaResult <- prcomp((freqs.ctr), scale=TRUE) 

loadings <- as.data.frame(pcaResult$rotation)

# now pull out extreme loadings
loads <- quantile(loadings$PC1, probs=c(0.005, 0.995))

# then match names of extreme loadings to pull out allele freqs from full data set
outliers <- loadings[which(loadings$PC1 < loads[1] | loadings$PC1 > loads[2]),]

row.names(af) <- paste(dat3$Chrom, dat3$Position, sep=":")

af.out <- af[(row.names(af) %in% row.names(outliers)),]

# and plot these allele freqs
pop.id <- unique(substr(colnames(af.out), 1,7))
mean.out <- as.data.frame(matrix(nrow=nrow(af.out),ncol=length(pop.id)))
colnames(mean.out) <- pop.id

for(i in 1:length(pop.id)){

    tmp_pop <- af.out[,grep(pop.id[i], colnames(af.out))]

    mean.out[,i] <- apply(tmp_pop, 1, mean)
}

png("~/reciprocal_t/figures/trt_allele_freq.png", height = 7, width = 9, res=300, units="in")

par(mfrow=c(2,3))
x <- mean.out$AAAA_F1 - mean.out$HHHH_F1
y <- mean.out$AAHH_F1 - mean.out$HHHH_F1
plot(x =x,
    y = y, 
    col=alpha("black", 0.1), bg=alpha("black", 0.05), pch=21,
    ylim=c(-1,1), xlim=c(-1, 1),
    xlab="AAAA - HHHH",
    ylab="(HH in AA) - HHHH",
    main="F1")
abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="red", lwd=2)
text(x=-0.7, y=0.7, paste("Corr = ", 
  round(cor.test(x,y, method="pearson")$estimate, 3)))

x <- mean.out$AAAA_F2 - mean.out$HHHH_F2
y <- mean.out$AAHH_F2 - mean.out$HHHH_F2
plot(x = x,
    y = y, 
    col=alpha("black", 0.1), bg=alpha("black", 0.05), pch=21,
    ylim=c(-1,1), xlim=c(-1, 1),
    xlab="AAAA - HHHH",
    ylab="(HH in AA) - HHHH",
    main="F2")
abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="red", lwd=2)
text(x=-0.7, y=0.7, paste("Corr = ", 
  round(cor.test(x,y, method="pearson")$estimate, 3)))

x <- mean.out$AAAA_F3 - mean.out$HHHH_F3
y <- mean.out$AAHH_F3 - mean.out$HHHH_F3
plot(x = x,
    y = y, 
    col=alpha("black", 0.1), bg=alpha("black", 0.05), pch=21,
    ylim=c(-1,1), xlim=c(-1, 1),
    xlab="AAAA - HHHH",
    ylab="(HH in AA) - HHHH",
    main="F3")
abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="red", lwd=2)
text(x=-0.7, y=0.7, paste("Corr = ", 
  round(cor.test(x,y, method="pearson")$estimate, 3)))

# and for AA in HH
x <- mean.out$HHHH_F1 - mean.out$AAAA_F1
y <- mean.out$HHAA_F1 - mean.out$AAAA_F1
plot(x =x,
    y = y, 
    col=alpha("black", 0.1), bg=alpha("black", 0.05), pch=21,
    ylim=c(-1,1), xlim=c(-1, 1),
    xlab="HHHH - AAAA",
    ylab="(AA in HH) - AAAA",
    main="F1")
abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="red", lwd=2)
text(x=-0.7, y=0.7, paste("Corr = ", 
  round(cor.test(x,y, method="pearson")$estimate, 3)))

x <- mean.out$HHHH_F2 - mean.out$AAAA_F2
y <- mean.out$HHAA_F2 - mean.out$AAAA_F2
plot(x = x,
    y = y, 
    col=alpha("black", 0.1), bg=alpha("black", 0.05), pch=21,
    ylim=c(-1,1), xlim=c(-1, 1),
    xlab="HHHH - AAAA",
    ylab="(AA in HH) - AAAA",
    main="F2")
abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="red", lwd=2)
text(x=-0.7, y=0.7, paste("Corr = ", 
  round(cor.test(x,y, method="pearson")$estimate, 3)))

x <- mean.out$HHHH_F3 - mean.out$AAAA_F3
y <- mean.out$HHAA_F3 - mean.out$AAAA_F3
plot(x = x,
    y = y, 
    col=alpha("black", 0.1), bg=alpha("black", 0.05), pch=21,
    ylim=c(-1,1), xlim=c(-1, 1),
    xlab="HHHH - AAAA",
    ylab="(AA in HH) - AAAA",
    main="F3")
abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="red", lwd=2)

text(x=-0.7, y=0.7, paste("Corr = ", 
  round(cor.test(x,y, method="pearson")$estimate, 3)))

dev.off()

#####
# decay in line correlation
#####
png("~/reciprocal_t/figures/line_allele_freq.png", height = 7, width = 9, res=300, units="in")

par(mfrow=c(2,3))
x <- mean.out$HHHH_F1 - mean.out$AAAA_F1
y <- mean.out$AAHH_F1 - mean.out$AAAA_F1
plot(x =x,
    y = y, 
    col=alpha("black", 0.1), bg=alpha("black", 0.05), pch=21,
    ylim=c(-1,1), xlim=c(-1, 1),
    xlab="HHHH - AAAA",
    ylab="(HH in AA) - AAAA",
    main="F1")
abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="red", lwd=2)
text(x=-0.7, y=0.7, paste("Corr = ", 
  round(cor.test(x,y, method="pearson")$estimate, 3)))

x <- mean.out$HHHH_F2 - mean.out$AAAA_F2
y <- mean.out$AAHH_F2 - mean.out$AAAA_F2
plot(x = x,
    y = y, 
    col=alpha("black", 0.1), bg=alpha("black", 0.05), pch=21,
    ylim=c(-1,1), xlim=c(-1, 1),
    xlab="HHHH - AAAA",
    ylab="(HH in AA) - AAAA",
    main="F2")
abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="red", lwd=2)
text(x=-0.7, y=0.7, paste("Corr = ", 
  round(cor.test(x,y, method="pearson")$estimate, 3)))

x <- mean.out$HHHH_F3 - mean.out$AAAA_F3
y <- mean.out$AAHH_F3 - mean.out$AAAA_F3
plot(x = x,
    y = y, 
    col=alpha("black", 0.1), bg=alpha("black", 0.05), pch=21,
    ylim=c(-1,1), xlim=c(-1, 1),
    xlab="HHHH - AAAA",
    ylab="(HH in AA) - AAAA",
    main="F3")
abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="red", lwd=2)
text(x=-0.7, y=0.7, paste("Corr = ", 
  round(cor.test(x,y, method="pearson")$estimate, 3)))

# and for AA in HH
x <- mean.out$AAAA_F1 - mean.out$HHHH_F1
y <- mean.out$HHAA_F1 - mean.out$HHHH_F1
plot(x =x,
    y = y, 
    col=alpha("black", 0.1), bg=alpha("black", 0.05), pch=21,
    ylim=c(-1,1), xlim=c(-1, 1),
    xlab="AAAA - HHHH",
    ylab="(AA in HH) - HHHH",
    main="F1")
abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="red", lwd=2)
text(x=-0.7, y=0.7, paste("Corr = ", 
  round(cor.test(x,y, method="pearson")$estimate, 3)))

x <- mean.out$AAAA_F2 - mean.out$HHHH_F2
y <- mean.out$HHAA_F2 - mean.out$HHHH_F2
plot(x = x,
    y = y, 
    col=alpha("black", 0.1), bg=alpha("black", 0.05), pch=21,
    ylim=c(-1,1), xlim=c(-1, 1),
    xlab="AAAA - HHHH",
    ylab="(AA in HH) - HHHH",
    main="F2")
abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="red", lwd=2)
text(x=-0.7, y=0.7, paste("Corr = ", 
  round(cor.test(x,y, method="pearson")$estimate, 3)))

x <- mean.out$AAAA_F3 - mean.out$HHHH_F3
y <- mean.out$HHAA_F3 - mean.out$HHHH_F3
plot(x = x,
    y = y, 
    col=alpha("black", 0.1), bg=alpha("black", 0.05), pch=21,
    ylim=c(-1,1), xlim=c(-1, 1),
    xlab="AAAA - HHHH",
    ylab="(AA in HH) - HHHH",
    main="F3")
abline(0,1, col="gray24", lty=2, lwd=2)
abline(h=0, col="gray24", lty=1, lwd=1)
abline(v=0, col="gray24", lty=1, lwd=1)
abline(lm(y~x), col="red", lwd=2)

text(x=-0.7, y=0.7, paste("Corr = ", 
  round(cor.test(x,y, method="pearson")$estimate, 3)))



dev.off()


