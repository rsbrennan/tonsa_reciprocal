# FIgure 3, popoolation_pi.R

library(stringr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(dplyr)

#locate the directory containing the files.
dir <- "~/reciprocal_t/analysis/popoolation"
#list.files(dir)
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



################
## plot mean het for each group. sep by treatments, etc.
################

het_mean <-as.data.frame(matrix(ncol=6, nrow=length(pops)))
colnames(het_mean) <- c("ID", "Generation", "comb", "Line", "Treatment", "Heterozygosity")
het_mean$ID <- pops
het_mean$comb <- substr(pops, 1,4)
het_mean$Treatment <- substr(pops, 1,2)
het_mean$Line <- substr(pops, 3,4)
het_mean$Generation <- substr(pops, 6,7)

for(i in 1:length(pops)){
    a <- median(d[[grep(pops[i], names(d))]]$pi, na.rm=TRUE)
    het_mean$Heterozygosity[which(het_mean$ID == pops[i])] <- a

}

het_mean$comb <- factor(het_mean$comb, levels = c("AAAA","HHAA", "HHHH", "AAHH"))

# anova

fit <- aov(Heterozygosity ~ comb*Generation, data=het_mean)
summary(fit)
stat <- TukeyHSD(fit, "comb:Generation")$`comb:Generation`
stat[which(stat[,4] < 0.1),]


totals <- het_mean %>%
    group_by(comb,Generation) %>%
             slice(which.max(Heterozygosity))

pb <- ggplot(data=het_mean, aes(x=comb, y=Heterozygosity, fill=Generation, shape=Generation)) +
  geom_point(position=position_dodge(width=0.95),aes(group=Generation), size=5.5, alpha=0.9) +

  scale_fill_manual(values=alpha(c("white","gray48", "black"), 0.8)) +
  scale_shape_manual(values=c(21, 21, 21)) +
  #geom_point(data=all_mean,aes(x=comb, y=Heterozygosity))+
  #geom_boxplot(aes(fill=Generation))+
  stat_summary(geom = "boxplot",
             fun.data = function(x) setNames(quantile(x, c(0.00, 0.25, 0.5, 0.75, 1)), c("ymin", "lower", "middle", "upper", "ymax")),
             position = "dodge" , show.legend=FALSE, aes(alpha=0.5))+
theme_base() +
scale_x_discrete(labels=c("AAAA" = "AM in AM",
                                "HHAA" = "AM in GH",
                                "HHHH" = "GH in GH",
                                "AAHH" = "GH in AM"))+
   xlab("")+ ylab(expression(paste("Median  ",pi))) +
   ylim(0.008,  0.0113) +
 #  stat_summary(geom = 'text', label = letters[1:12], fun.y = (max), position = position_dodge(width=0.95))
    geom_text(aes(x=comb, y=Heterozygosity+0.0002,
        label =c("AB", "AB", "AB", "AB", "A", "A", "AB", "ABC", "BC","AB", "AB", "C"),
        fill = NULL), data = totals,
        position = position_dodge(width=0.95), size=4)+
     theme(axis.text.x = element_text(angle=45, hjust=1, size=16),
        axis.title.y = element_text( size=16))

pb

pdf("~/reciprocal_t/figures/Fig_4_pi.pdf", height=5, width=7)

pb

dev.off()
