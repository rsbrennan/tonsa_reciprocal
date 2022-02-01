
# pi and selection. after grc conference.

library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(dplyr)
library(Rmisc)
library(reshape2)
library(ggthemes)
library(ggpubr)

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

for(i in 1:length(d)){
 d[[i]]$gp <- names(d)[i]

}

pops <- names(d)

###############################################
########
######## plot het for each group
########
################################################

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


pb <- ggplot(data=het_mean, aes(x=comb, y=Heterozygosity, fill=Generation, shape=Generation)) +
  geom_point(position=position_dodge(width=0.95),aes(group=Generation), size=5.5, alpha=0.9) +

  scale_fill_manual(values=alpha(c("white","gray48", "black"), 0.8)) +
  scale_shape_manual(values=c(21, 21, 21)) +
  #geom_point(data=all_mean,aes(x=comb, y=Heterozygosity))+
  #geom_boxplot(aes(fill=Generation))+
  stat_summary(geom = "boxplot",
             fun.data = function(x) setNames(quantile(x, c(0.00, 0.25, 0.5, 0.75, 1)), c("ymin", "lower", "middle", "upper", "ymax")),
             position = "dodge" , show.legend=FALSE, aes(alpha=0.5))+
theme_classic() +
scale_x_discrete(labels=c("AAAA" = "AM in AM",
                                "HHAA" = "AM in GH",
                                "HHHH" = "GH in GH",
                                "AAHH" = "GH in AM"))+
   xlab("")+ ylab(expression(paste("Median  ",pi))) +
   ylim(0.008,  0.0113) +
 #  stat_summary(geom = 'text', label = letters[1:12], fun.y = (max), position = position_dodge(width=0.95))
    geom_text(aes(x=comb, y=Heterozygosity+0.0002,
        label =c("AB", "AB", "AB", "AB", "A", "A", "AB", "ABC", "BC","AB", "AB", "C"),
        fill = NULL), data = het_mean,
        position = position_dodge(width=0.95), size=4)+
     theme(axis.text.x = element_text(angle=45, hjust=1, size=16),
        axis.title.y = element_text( size=14))


######################
#
# where is diversity being lost?
#
######################

# write output to look for overlap. will use bedtools.

dpi <- bind_rows(d, .id = "column_label")

dout <- dpi[(which(dpi$pi != "na")),]

dat <- read.csv("~/reciprocal_t/analysis/AF_change.txt", sep="\t")

# split snp into chr and snp
dat$CHR <- str_split_fixed(dat$SNP, ":", 2)[,1]
dat$POS <- str_split_fixed(dat$SNP, ":", 2)[,2]
write.table(cbind(dat$CHR, as.numeric(dat$POS)-1, dat$POS, dat$sig_all),
        file="~/reciprocal_t/analysis/sig.bed",
        quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

dout$start <- dout$position-50
dout$stop <- dout$position+51

# add pi measure to this. and group
write.table(cbind(dout$gene, dout$start, dout$stop, dout[4:8]),
        file="~/reciprocal_t/analysis/pi.bed",
        quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

# the following is run in bash

#sort -k1,1 -k2,2n pi.bed > pi.sorted.bed
#sort -k1,1 -k2,2n sig.bed > sig.sorted.bed
#
#bedtools intersect -wa -wb \
#    -a pi.sorted.bed \
#    -b sig.sorted.bed \
#    -sorted > ~/reciprocal_t/analysis/pi.overlp.bed


dat <- read.csv("~/reciprocal_t/analysis/pi.overlp.bed", header=FALSE, sep="\t")

colnames(dat) <- c("CHR", "start", "stop", "cov", "nsnp", "pi",
                    "pi_snp", "gp", "gene", "snp_start", "snp_stop", "sig")
dat$trt <- str_split_fixed(dat$gp, "_", 2)[,2]
dat$trt <- substr(dat$gp, 1,4)

################
###
### delta pi. to look for enrichment
###
################

# use list of dfs d

# merge the data sets.
# the look for regions with drops.

HHAA_r1 <- merge(d[["AAAA_F1_REP1"]], d[["HHAA_F3_REP1"]], by="snp")
HHAA_r2 <- merge(d[["AAAA_F1_REP2"]], d[["HHAA_F3_REP2"]], by="snp")
HHAA_r3 <- merge(d[["AAAA_F1_REP3"]], d[["HHAA_F3_REP3"]], by="snp")
HHAA_r4 <- merge(d[["AAAA_F1_REP4"]], d[["HHAA_F3_REP4"]], by="snp")

AAHH_r1 <- merge(d[["HHHH_F1_REP1"]], d[["AAHH_F3_REP1"]], by="snp")
AAHH_r2 <- merge(d[["HHHH_F1_REP2"]], d[["AAHH_F3_REP2"]], by="snp")
AAHH_r3 <- merge(d[["HHHH_F1_REP3"]], d[["AAHH_F3_REP3"]], by="snp")
AAHH_r4 <- merge(d[["HHHH_F1_REP4"]], d[["AAHH_F3_REP4"]], by="snp")

HHAA_r1$pi.HHAA_r1 <- HHAA_r1$pi.y-HHAA_r1$pi.x
HHAA_r2$pi.HHAA_r2 <- HHAA_r2$pi.y-HHAA_r2$pi.x
HHAA_r3$pi.HHAA_r3 <- HHAA_r3$pi.y-HHAA_r3$pi.x
HHAA_r4$pi.HHAA_r4 <- HHAA_r4$pi.y-HHAA_r4$pi.x

AAHH_r1$pi.AAHH_r1 <- AAHH_r1$pi.y-AAHH_r1$pi.x
AAHH_r2$pi.AAHH_r2 <- AAHH_r2$pi.y-AAHH_r2$pi.x
AAHH_r3$pi.AAHH_r3 <- AAHH_r3$pi.y-AAHH_r3$pi.x
AAHH_r4$pi.AAHH_r4 <- AAHH_r4$pi.y-AAHH_r4$pi.x

HHAA_r1 <- HHAA_r1[!is.na(HHAA_r1$pi.HHAA_r1),]
HHAA_r2 <- HHAA_r2[!is.na(HHAA_r2$pi.HHAA_r2),]
HHAA_r3 <- HHAA_r3[!is.na(HHAA_r3$pi.HHAA_r3),]
HHAA_r4 <- HHAA_r4[!is.na(HHAA_r4$pi.HHAA_r4),]

AAHH_r1 <- AAHH_r1[!is.na(AAHH_r1$pi.AAHH_r1),]
AAHH_r2 <- AAHH_r2[!is.na(AAHH_r2$pi.AAHH_r2),]
AAHH_r3 <- AAHH_r3[!is.na(AAHH_r3$pi.AAHH_r3),]
AAHH_r4 <- AAHH_r4[!is.na(AAHH_r4$pi.AAHH_r4),]

# there are dup cols in these merges, but doesn't matter. ignore. not worth fixing
all.HHAA <- merge(merge(merge(HHAA_r1, HHAA_r2, by="snp"),HHAA_r3, by="snp"),HHAA_r4, by="snp")
all.AAHH <- merge(merge(merge(AAHH_r1,AAHH_r2, by="snp"),AAHH_r3, by="snp"),AAHH_r4, by="snp")

#all.m <- merge(merge(merge(all.HHAA, all.AAHH, by="snp", all=T),all.AAAA, by="snp", all=T), all.HHHH, by="snp", all=T)
all.m <- merge(all.HHAA, all.AAHH, by="snp")

nrow(all.m)
# [1] 9177


################
# save output for GO enrichment.
################

new.dat <- data.frame(
            snp=all.m$snp,
            gene=all.m$gene.x.x.x,
            AAHH_r1 = all.m$pi.AAHH_r1,
            AAHH_r2 = all.m$pi.AAHH_r2,
            AAHH_r3 = all.m$pi.AAHH_r3,
            AAHH_r4 = all.m$pi.AAHH_r4,
            HHAA_r1 = all.m$pi.HHAA_r1,
            HHAA_r2 = all.m$pi.HHAA_r2,
            HHAA_r3 = all.m$pi.HHAA_r3,
            HHAA_r4 = all.m$pi.HHAA_r4)

df.melt <- melt(new.dat, id.vars = c("snp", "gene"))
df.melt$gp <- substr(df.melt$variable, 1,4)

data <-
  df.melt %>%
  group_by(gp, gene) %>%
  dplyr::summarise(mean_pi = mean(value))

aahh <- data[which(data$gp == "AAHH"),]
hhaa <- data[which(data$gp == "HHAA"),]


write.table(file="~/reciprocal_t/analysis/GO_enrich/pi_aahh.txt", aahh[,c(2,3)], col.names=TRUE,
    row.names=FALSE, quote=FALSE,sep=",")
write.table(file="~/reciprocal_t/analysis/GO_enrich/pi_hhaa.txt", hhaa[,c(2,3)], col.names=TRUE,
    row.names=FALSE, quote=FALSE,sep=",")


# window level
df.melt <- melt(new.dat, id.vars = c("snp", "gene"))
df.melt$gp <- substr(df.melt$variable, 1,4)

deltapi <-
  df.melt %>%
  group_by(gp, snp, gene) %>%
  dplyr::summarise(mean_pi = mean(value, ignore.na=TRUE))

# remember that dat contains the overlaps.
    # so could merge deltapi with dat to see significance

overlap <- (merge(x=deltapi, y=dat, by.x="snp", by.y="pi_snp", all.x=TRUE))
# messy, columns that are carried over. remove them
clean <- overlap[ , c("snp","gp.x", "gene.x", "mean_pi", "CHR","start", "stop","sig", "snp_start", "snp_stop")]

colnames(clean) <- c("snp","gp", "gene", "mean_pi", "CHR","start", "stop","sig", "snp_start", "snp_stop")

# and if no snp, this means sig=false, bc didn't shift.
clean$sig[is.na(clean$sig)] <- FALSE
clean$sig[is.na(clean$sig)] <- FALSE

deduped.data <- clean %>% distinct

deduped.data$gp <- factor(deduped.data$gp, levels = c("HHAA", "AAHH"))

c <- ggplot(deduped.data, aes(x=gp, y=mean_pi, color=sig, group=sig)) +
     stat_summary(fun.data = mean_cl_normal,position=position_dodge(0.5), geom = "errorbar") +
     stat_summary(fun.data = mean_cl_normal,position=position_dodge(0.5), geom = "point", size=4) +
     theme_classic() +
     scale_colour_manual(values = c("gray48", "black"), name = " ", labels = c("Non-adaptive", "Adaptive")) +
     theme(axis.text.x = element_text(angle=45, hjust=1, size=16),
          axis.title.y = element_text( size=16)) +
     scale_x_discrete(labels=c( "HHAA" = "AM in GH",
                                "AAHH" = "GH in AM")) +
     xlab("")+
     ylab(expression(paste("Mean change in ",pi))) +
     theme(axis.title.y=element_text(size=14))
    # ylab(expression(paste("Mean loss of  ",pi, "\n(100 bp windows")))

ggsave("~/reciprocal_t/figures/delta_pi.png",
        plot =ggarrange(pb, c, widths = c(2,1.3), ncol=2, nrow=1, labels="AUTO"), 
        width=10, height=4)

ggsave("~/reciprocal_t/figures/delta_pi.pdf",
    plot =ggarrange(pb, c, widths = c(2,1.3), ncol=2, nrow=1, labels="AUTO"),
        width=10, height=4)




win.aahh <- deduped.data[which(deduped.data$gp == "AAHH"),]

wilcox.test(win.aahh$mean_pi[which(win.aahh$sig == TRUE)],
              win.aahh$mean_pi[which(win.aahh$sig == FALSE)],
              conf.int = TRUE, alternative= "less")


win.hhaa <- deduped.data[which(deduped.data$gp == "HHAA"),]

wilcox.test(win.hhaa$mean_pi[which(win.hhaa$sig == TRUE)],
              win.hhaa$mean_pi[which(win.hhaa$sig == FALSE)],
              conf.int = TRUE, alternative= "less")

p.adjust(c(0.002297, 0.0534), method="bonferroni")

############################################################################################################
############################################################################################################
# do AMam VS GHgh show changes in pi?
############################################################################################################
############################################################################################################


################
###
### delta pi. to look for enrichment
###
################

# use list of dfs d

# merge the data sets.
# the look for regions with drops.

HHHH_r1 <- merge(d[["AAAA_F1_REP1"]], d[["HHHH_F1_REP1"]], by="snp")
HHHH_r2 <- merge(d[["AAAA_F1_REP2"]], d[["HHHH_F1_REP2"]], by="snp")
HHHH_r3 <- merge(d[["AAAA_F1_REP3"]], d[["HHHH_F1_REP3"]], by="snp")
HHHH_r4 <- merge(d[["AAAA_F1_REP4"]], d[["HHHH_F1_REP4"]], by="snp")

HHHH_r1$pi.HHHH_r1 <- HHHH_r1$pi.y-HHHH_r1$pi.x
HHHH_r2$pi.HHHH_r2 <- HHHH_r2$pi.y-HHHH_r2$pi.x
HHHH_r3$pi.HHHH_r3 <- HHHH_r3$pi.y-HHHH_r3$pi.x
HHHH_r4$pi.HHHH_r4 <- HHHH_r4$pi.y-HHHH_r4$pi.x

HHHH_r1 <- HHHH_r1[!is.na(HHHH_r1$pi.HHHH_r1),]
HHHH_r2 <- HHHH_r2[!is.na(HHHH_r2$pi.HHHH_r2),]
HHHH_r3 <- HHHH_r3[!is.na(HHHH_r3$pi.HHHH_r3),]
HHHH_r4 <- HHHH_r4[!is.na(HHHH_r4$pi.HHHH_r4),]

# there are dup cols in these merges, but doesn't matter. ignore. not worth fixing
all.HHHH <- merge(merge(merge(HHHH_r1, HHHH_r2, by="snp"),HHHH_r3, by="snp"),HHHH_r4, by="snp")

#all.m <- merge(merge(merge(all.HHAA, all.AAHH, by="snp", all=T),all.AAAA, by="snp", all=T), all.HHHH, by="snp", all=T)

nrow(all.HHHH)
# [1] 11648

all.m <- merge(all.HHAA, all.AAHH, by="snp")
#all.m <- merge(merge(merge(all.HHAA, all.AAHH, by="snp", all=T),all.AAAA, by="snp", all=T), all.HHHH, by="snp", all=T)
all.combo <- merge(all.m, all.HHHH, by="snp")

nrow(all.combo)
# 9177

################
# save output for GO enrichment.
################

new.dat <- data.frame(
      snp=all.combo$snp,
      gene=all.combo$gene.x.x.x,
      HHHH_r1 = all.combo$pi.HHHH_r1,
      HHHH_r2 = all.combo$pi.HHHH_r2,
      HHHH_r3 = all.combo$pi.HHHH_r3,
      HHHH_r4 = all.combo$pi.HHHH_r4,
      AAHH_r1 = all.combo$pi.AAHH_r1,
      AAHH_r2 = all.combo$pi.AAHH_r2,
      AAHH_r3 = all.combo$pi.AAHH_r3,
      AAHH_r4 = all.combo$pi.AAHH_r4,
      HHAA_r1 = all.combo$pi.HHAA_r1,
      HHAA_r2 = all.combo$pi.HHAA_r2,
      HHAA_r3 = all.combo$pi.HHAA_r3,
      HHAA_r4 = all.combo$pi.HHAA_r4)

df.melt <- melt(new.dat, id.vars = c("snp", "gene"))
df.melt$gp <- substr(df.melt$variable, 1,4)

data <-
  df.melt %>%
  group_by(gp, gene) %>%
  dplyr::summarise(mean_pi = mean(value))

hhhh <- data[which(data$gp == "HHHH"),]

write.table(file="~/reciprocal_t/analysis/GO_enrich/pi_hhhh.txt", hhhh[,c(2,3)], col.names=TRUE,
    row.names=FALSE, quote=FALSE,sep=",")


# window level
df.melt <- melt(new.dat, id.vars = c("snp", "gene"))
df.melt$gp <- substr(df.melt$variable, 1,4)

deltapi <-
  df.melt %>%
  group_by(gp, snp, gene) %>%
  dplyr::summarise(mean_pi = mean(value, ignore.na=TRUE))

# remember that dat contains the overlaps.
  # so could merge deltapi with dat to see significance

overlap <- (merge(x=deltapi, y=dat, by.x="snp", by.y="pi_snp", all.x=TRUE))
# messy, columns that are carried over. remove them
clean <- overlap[ , c("snp","gp.x", "gene.x", "mean_pi", "CHR","start", "stop","sig", "snp_start", "snp_stop")]

colnames(clean) <- c("snp","gp", "gene", "mean_pi", "CHR","start", "stop","sig", "snp_start", "snp_stop")

# and if no snp, this means sig=false, bc didn't shift.
clean$sig[is.na(clean$sig)] <- FALSE
clean$sig[is.na(clean$sig)] <- FALSE

deduped.data <- clean %>% distinct

deduped.data$gp <- factor(deduped.data$gp, levels = c("HHAA","HHHH", "AAHH"))

c <- ggplot(deduped.data, aes(x=gp, y=mean_pi, color=sig, group=sig)) +
     stat_summary(fun.data = mean_cl_normal,position=position_dodge(0.5), geom = "errorbar") +
     stat_summary(fun.data = mean_cl_normal,position=position_dodge(0.5), geom = "point", size=4) +
     theme_classic() +
     scale_colour_manual(values = c("gray48", "black"), name = " ", labels = c("Non-adaptive", "Adaptive")) +
     theme(axis.text.x = element_text(angle=45, hjust=1, size=16),
        axis.title.y = element_text( size=16)) +
     scale_x_discrete(labels=c( "HHHH" = "GH in GH",
                                "HHAA" = "AM in GH",
                                "AAHH" = "GH in AM")) +
     xlab("")+
     ylab(expression(paste("Mean change in ",pi))) +
     theme(axis.title.y=element_text(size=14)) 
         # ylab(expression(paste("Mean loss of  ",pi, "\n(100 bp windows")))

#ggsave("~/reciprocal_t/figures/delta_pi.png",
 #   plot =ggarrange(pb, c, widths = c(2,1.3), ncol=2, nrow=1, labels="AUTO"), 
 #       width=10, height=4)

write.table(deduped.data, file="~/reciprocal_t/analysis/Fig5_b.txt", row.names=F, quote=F)

ggsave("~/reciprocal_t/figures/delta_pi_gh.pdf",
    plot =ggarrange(pb,c, widths = c(2,1.3), ncol=2, nrow=1, labels="AUTO"),
        width=10, height=4)




win.hhhh <- deduped.data[which(deduped.data$gp == "HHHH"),]

wilcox.test(win.hhhh$mean_pi[which(win.hhhh$sig == TRUE)],
        win.hhhh$mean_pi[which(win.hhhh$sig == FALSE)],
        conf.int = TRUE, alternative= "less")
# 0.02177

win.hhaa <- deduped.data[which(deduped.data$gp == "HHAA"),]

wilcox.test(win.hhaa$mean_pi[which(win.hhaa$sig == TRUE)],
        win.hhaa$mean_pi[which(win.hhaa$sig == FALSE)],
        conf.int = TRUE, alternative= "less")
# 0.0534

win.aahh <- deduped.data[which(deduped.data$gp == "AAHH"),]

wilcox.test(win.aahh$mean_pi[which(win.aahh$sig == TRUE)],
        win.aahh$mean_pi[which(win.aahh$sig == FALSE)],
        conf.int = TRUE, alternative= "less")
# 0.002297

p.adjust(c(0.02177, 0.0534,0.002297), method = "bonferroni")
# [1] 0.065310 0.160200 0.006891
