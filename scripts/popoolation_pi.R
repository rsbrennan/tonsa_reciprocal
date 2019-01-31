
library(stringr)
library(tidyr)

#locate the directory containing the files.
dir <- "~/reciprocal_t/analysis/popoolation"
list.files(dir)
files <- file.path(dir, list.files(dir))

files <- files[grep("100bp.pi", files)]
files <- files[grep("params", files, invert=TRUE)]
d <- lapply(files, read.table)
names(d) <- (str_split_fixed(files, "/", 5)[,5] %>% str_split_fixed( "[.]", 3))[,1]

d <- lapply(d, function(x) {
  colnames(x) <- c("gene", "position", "n_variants", "prop_covered", "pi")
  x
})

d <- lapply(d, function(x) {
  x$pi <- as.numeric(as.character(x$pi))
  x
})

d <- lapply(d, function(x) {
  x$snp <- paste(x$gene, x$position, sep="_")
  x
})

# count number of windows for each sample:
lapply(d, function(x) {
 length(which(x$pi != "na"))
})

pops <- names(d)

# hist plots

par(mfrow = c(1, 1), mar=c(3, 3, 3, 1), mgp=c(3, 1, 0), las=0)

plot(1, type="n", ylim=c(0,50), xlim=c(0,0.3),
    lwd=2, xaxt="n",yaxt="n", main="", xlab="", ylab="")
axis(1, mgp=c(2, .5, 0), cex.axis=0.7) # second is tick mark labels
axis(2, mgp=c(2, .5, 0), cex.axis=0.7)

title(xlab="heterozygosity", line=2, cex.lab=0.9)
title(ylab="density", line=2, cex.lab=0.9)
title(main="heterozygosity: transcript avg", line=1, cex.lab=0.9)

for(i in 1:length(d)){
    if(substr(names(d)[[i]], 1,4) == "AAAA"){
        lines(density(d[[i]]$pi, na.rm=TRUE), lwd=2, col=alpha("black", 0.2))
    }
    if(substr(names(d)[[i]], 1,4) == "HHHH"){
        lines(density(d[[i]]$pi, na.rm=TRUE), lwd=2, col=alpha("red", 0.2))
    }
    if(substr(names(d)[[i]], 1,4) == "HHAA"){
        lines(density(d[[i]]$pi, na.rm=TRUE), lty=2, lwd=2, col=alpha("black", 0.2))
    }
    if(substr(names(d)[[i]], 1,4) == "AAHH"){
        lines(density(d[[i]]$pi, na.rm=TRUE), lty=2, lwd=2, col=alpha("red", 0.2))
    }
}


legend("topright", c("AAAA", "HHHH", "AA in HH", "HH in AA"),
    col=c("black", "red", "black", "red"), lwd=2, lty=c(1,1,2,2), cex=0.8 )



## plot mean het for each group. sep by treatments, etc. 

het_mean <-as.data.frame(matrix(ncol=6, nrow=length(pops)))
colnames(het_mean) <- c("ID", "Generation", "comb", "Line", "Treatment", "Heterozygosity")
het_mean$ID <- pops
het_mean$comb <- substr(pops, 1,4)
het_mean$Treatment <- substr(pops, 1,2)
het_mean$Line <- substr(pops, 3,4)
het_mean$Generation <- substr(pops, 6,7)

for(i in 1:length(pops)){
    a <- mean(d[[grep(pops[i], names(d))]]$pi, na.rm=TRUE)
    het_mean$Heterozygosity[which(het_mean$ID == pops[i])] <- a

}

het_mean$comb <- factor(het_mean$comb, levels = c("AAAA","HHAA", "HHHH", "AAHH"))

pa <- ggplot(data=het_mean, aes(x=comb, y=Heterozygosity, fill=Generation, shape=Generation)) +
  geom_point(position=position_dodge(width=0.85),aes(group=Generation), size=4, alpha=0.9) +

  scale_fill_manual(values=alpha(c("darkolivegreen3","gray60", "darkgoldenrod3"), 0.8)) +
  scale_shape_manual(values=c(21, 24, 22)) +
  #geom_point(data=all_mean,aes(x=comb, y=Heterozygosity))+
  #geom_boxplot(aes(fill=Generation))+
  stat_summary(geom = "boxplot", 
             fun.data = function(x) setNames(quantile(x, c(0.00, 0.25, 0.5, 0.75, 1)), c("ymin", "lower", "middle", "upper", "ymax")), 
             position = "dodge" , show.legend=FALSE)+
    theme_bw()+
   scale_x_discrete(labels=c("AAAA" = "AA in AA", 
                                "HHAA" = "AA in HH",
                                "HHHH" = "HH in HH",
                                "AAHH" = "HH in AA"))+
   xlab("")+ ggtitle("pi: mean")



for(i in 1:length(pops)){
    a <- median(d[[grep(pops[i], names(d))]]$pi, na.rm=TRUE)
    het_mean$Heterozygosity[which(het_mean$ID == pops[i])] <- a

}



# solution: add letters for post-hoc Tukey results above upper whiskers

library(dplyr)
totals <- het_mean %>%
    group_by(comb,Generation) %>%
             slice(which.max(Heterozygosity))

pb <- ggplot(data=het_mean, aes(x=comb, y=Heterozygosity, fill=Generation, shape=Generation)) +
  geom_point(position=position_dodge(width=0.95),aes(group=Generation), size=4, alpha=0.9) +

  scale_fill_manual(values=alpha(c("darkolivegreen3","gray60", "darkgoldenrod3"), 0.8)) +
  scale_shape_manual(values=c(21, 24, 22)) +
  #geom_point(data=all_mean,aes(x=comb, y=Heterozygosity))+
  #geom_boxplot(aes(fill=Generation))+
  stat_summary(geom = "boxplot", 
             fun.data = function(x) setNames(quantile(x, c(0.00, 0.25, 0.5, 0.75, 1)), c("ymin", "lower", "middle", "upper", "ymax")), 
             position = "dodge" , show.legend=FALSE)+
    theme_bw()+
   scale_x_discrete(labels=c("AAAA" = "AA in AA", 
                                "HHAA" = "AA in HH",
                                "HHHH" = "HH in HH",
                                "AAHH" = "HH in AA"))+
   xlab("")+ ylab(expression(paste("Median  ",pi))) +
   ylim(0.008,  0.0113) +
 #  stat_summary(geom = 'text', label = letters[1:12], fun.y = (max), position = position_dodge(width=0.95))
    geom_text(aes(x=comb, y=Heterozygosity+0.0002,
        label =c("AB", "AB", "AB", "AB", "A", "A", "AB", "ABC", "BC","AB", "AB", "C"),
        fill = NULL), data = totals,
        position = position_dodge(width=0.95), size=4)+
     theme(axis.text.x = element_text(angle=45, hjust=1, size=16),
        axis.title.y = element_text( size=16)) 


png("~/reciprocal_t/figures/pi.png", res=300, height=5, width=7, units="in")


pb

dev.off()


#library(ggpubr)
ggarrange(pa, pb, common.legend=TRUE)


# anova

fit <- aov(Heterozygosity ~ comb*Generation, data=het_mean)
summary(fit)
stat <- TukeyHSD(fit, "comb:Generation")$`comb:Generation`
stat[which(stat[,4] < 0.1),]


From mean:
                         diff          lwr           upr        p adj
AAHH:F3-AAAA:F1 -0.0008801273 -0.001739946 -2.030872e-05 0.0408959252
AAHH:F3-HHAA:F1 -0.0009794685 -0.001839287 -1.196500e-04 0.0145225165
AAHH:F3-HHHH:F1 -0.0008755688 -0.001735387 -1.575020e-05 0.0427981693
AAHH:F3-AAHH:F1 -0.0008053997 -0.001665218  5.441884e-05 0.0838741880
AAHH:F3-AAAA:F2 -0.0008356531 -0.001695472  2.416553e-05 0.0631638507
HHHH:F3-HHAA:F2 -0.0008499798 -0.001709798  9.838736e-06 0.0550318561
AAHH:F3-HHAA:F2 -0.0013404147 -0.002200233 -4.805961e-04 0.0002171632
AAHH:F3-AAHH:F2 -0.0010674734 -0.001927292 -2.076548e-04 0.0054765217
AAHH:F3-AAAA:F3 -0.0009465849 -0.001806404 -8.676635e-05 0.0206384095
AAHH:F3-HHAA:F3 -0.0012630088 -0.002122827 -4.031902e-04 0.0005535832

Median:
                         diff          lwr           upr        p adj
AAHH:F3-AAAA:F1 -0.0010340375 -0.002022648 -4.542737e-05 0.0336965744
AAHH:F3-HHAA:F1 -0.0012388518 -0.002227462 -2.502416e-04 0.0048893903
AAHH:F3-HHHH:F1 -0.0009475979 -0.001936208  4.101226e-05 0.0704417695
AAHH:F3-AAHH:F1 -0.0010329975 -0.002021608 -4.438737e-05 0.0340079284
AAHH:F3-AAAA:F2 -0.0010524880 -0.002041098 -6.387787e-05 0.0285891071
HHHH:F3-HHAA:F2 -0.0010737177 -0.002062328 -8.510762e-05 0.0235986536
AAHH:F3-HHAA:F2 -0.0015070006 -0.002495611 -5.183905e-04 0.0003114930
AAHH:F3-AAHH:F2 -0.0011603133 -0.002148923 -1.717031e-04 0.0105064761
AAHH:F3-AAAA:F3 -0.0010232169 -0.002011827 -3.460674e-05 0.0370670656
HHHH:F3-HHAA:F3 -0.0011032079 -0.002091818 -1.145977e-04 0.0179969699
AAHH:F3-HHAA:F3 -0.0015364907 -0.002525101 -5.478806e-04 0.0002282283





for(i in 1:length(pops)){
    d[[i]]$id <- names(d)[i]
}


library(dplyr)
mdat <- bind_rows(d, .id = "id")

mdat$Line <- substr(mdat$id, 3,4)
mdat$Treatment <- substr(mdat$id, 1,2)
mdat$Rep <- substr(mdat$id,9,12)
mdat$Gen <- substr(mdat$id,6,7)
mdat$Group <- substr(mdat$id,1,4)


fit <- aov(pi ~ Group*Gen, data=mdat)
summary(fit)
stat <- TukeyHSD(fit, "Group:Gen")$`Group:Gen`
stat[which(stat[,4] < 0.01),]

library(nlme)
library(lme4)
mdat$Group <- as.factor(mdat$Group)
mdat$Gen <- as.factor(mdat$Gen)
mdat$Rep <- as.factor(mdat$Rep)
mdat2 <- mdat[!is.na(mdat$pi),]

m1 <- lmer(pi ~ Group*Gen+(1|Rep), data=mdat2,na.action=na.exclude, REML=FALSE)
m2 <- lmer(pi ~ Group+(1|Rep), data=mdat2,na.action=na.exclude, REML=FALSE)
m3 <- lmer(pi ~ 1+(1|Rep), data=mdat2,na.action=na.exclude, REML=FALSE)
summary(m1)

anova(m2,m1)
anova(m3,m2)

