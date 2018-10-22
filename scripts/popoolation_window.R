library(stringr)
library(tidyr)

#locate the directory containing the files. 
dir <- "~/reciprocal_t/analysis"
list.files(dir)
files <- file.path(dir, list.files(dir))

files <- files[grep("1kb.pi", files)]
files <- files[grep("params", files, invert=TRUE)]
d <- lapply(files, read.table)
names(d) <- (str_split_fixed(files, "/", 4)[,4] %>% str_split_fixed( "[.]", 3))[,1]

d <- lapply(d, function(x) {
  colnames(x) <- c("gene", "position", "n_variants", "prop_covered", "pi")
  x      
})

d <- lapply(d, function(x) {
  x$pi <- as.numeric(as.character(x$pi))
  x      
})

d <- lapply(d, function(x) {
  x$id <- paste(x$gene, x$position, sep="_")
  x      
})

# count number of windows for each sample:
lapply(d, function(x) {
 length(which(x$pi != "na"))      
})

##### mean

mean(d[['AAAA_F1_REP1']]$pi, na.rm=TRUE)
mean(d[['AAAA_F1_REP2']]$pi, na.rm=TRUE)
mean(d[['AAAA_F1_REP3']]$pi, na.rm=TRUE)
mean(d[['AAAA_F1_REP4']]$pi, na.rm=TRUE)
mean(d[['HHHH_F1_REP1']]$pi, na.rm=TRUE)
mean(d[['HHHH_F1_REP2']]$pi, na.rm=TRUE)
mean(d[['HHHH_F1_REP3']]$pi, na.rm=TRUE)
mean(d[['HHHH_F1_REP4']]$pi, na.rm=TRUE)

#######

png("~/reciprocal_t/figures/1kb.pi.hist.png", res=300, height=6, width=6, units="in")

par(mfrow = c(1, 1), mar=c(3, 3, 3, 1), mgp=c(3, 1, 0), las=0)

plot(density( d[['AAAA_F1_REP1']]$pi, na.rm=TRUE), 
    lwd=2, xaxt="n",yaxt="n", main="")
lines(density(d[['AAAA_F1_REP2']]$pi, na.rm=TRUE), lwd=2)
lines(density(d[['AAAA_F1_REP3']]$pi, na.rm=TRUE), lwd=2)
lines(density(d[['AAAA_F1_REP4']]$pi, na.rm=TRUE), lwd=2)
lines(density(d[['HHHH_F1_REP1']]$pi, na.rm=TRUE), lwd=2, col="red")
lines(density(d[['HHHH_F1_REP2']]$pi, na.rm=TRUE), lwd=2, col="red")
lines(density(d[['HHHH_F1_REP3']]$pi, na.rm=TRUE), lwd=2, col="red")
lines(density(d[['HHHH_F1_REP4']]$pi, na.rm=TRUE), lwd=2, col="red")
legend("topright", c("AAAA", "HHHH"), 
    col=c("black", "red"), lwd=2, cex=0.8 )

axis(1, mgp=c(2, .5, 0), cex.axis=0.7) # second is tick mark labels
axis(2, mgp=c(2, .5, 0), cex.axis=0.7)

title(xlab="pi", line=2, cex.lab=0.9)
title(ylab="density", line=2, cex.lab=0.9)
title(main="pi: 1kb windows", line=1, cex.lab=0.9)

dev.off()

a1 <- quantile(d[['AAAA_F1_REP1']]$pi, na.rm=TRUE, c(0.01, 0.05))
a2 <- quantile(d[['AAAA_F1_REP2']]$pi, na.rm=TRUE, c(0.01, 0.05))
a3 <- quantile(d[['AAAA_F1_REP3']]$pi, na.rm=TRUE, c(0.01, 0.05))
a4 <- quantile(d[['AAAA_F1_REP4']]$pi, na.rm=TRUE, c(0.01, 0.05))
h1 <- quantile(d[['HHHH_F1_REP1']]$pi, na.rm=TRUE, c(0.01, 0.05))
h2 <- quantile(d[['HHHH_F1_REP2']]$pi, na.rm=TRUE, c(0.01, 0.05))
h3 <- quantile(d[['HHHH_F1_REP3']]$pi, na.rm=TRUE, c(0.01, 0.05))
h4 <- quantile(d[['HHHH_F1_REP4']]$pi, na.rm=TRUE, c(0.01, 0.05))

gdat_05 <- data.frame(line=c(rep("AAAA", 4), rep("HHHH", 4)),
            pi=c(a1[2], a2[2], a3[2], a4[2], h1[2], h2[2], h3[2], h4[2]))

gdat <- data.frame(line=c(rep("AAAA", 4), rep("HHHH", 4)),
            pi=c(a1[1], a2[1], a3[1], a4[1], h1[1], h2[1], h3[1], h4[1]))

library(ggplot2)

png("~/reciprocal_t/figures/1kb.pi01_quantile.png",width = 6, height = 6, res=300, units="in")

ggplot(gdat, aes(x=line, y=pi, color=line)) +
    theme_bw() +
    stat_summary(fun.data="mean_se", mapping=aes(group=line), size=2.5,alpha = 0.5, show.legend=FALSE)+
    geom_jitter(width=0.1, size=3.5) +
    scale_color_manual(values=c("black", "firebrick3")) +
    ggtitle("1kb pi: 0.01 quantile") +
    theme(axis.text.x=element_text(size=14),
        axis.title.x=element_blank()) +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5))

dev.off()

png("~/reciprocal_t/figures/1kb.pi05_quantile.png",width = 6, height = 6, res=300, units="in")
ggplot(gdat_05, aes(x=line, y=pi, color=line)) +
    theme_bw() +
    stat_summary(fun.data="mean_se", mapping=aes(group=line), size=2.5,alpha = 0.5, show.legend=FALSE)+
    geom_jitter(width=0.1, size=3.5) +
    scale_color_manual(values=c("black", "firebrick3")) +
    ggtitle("1kb pi: 0.05 quantile") +
    theme(axis.text.x=element_text(size=14),
        axis.title.x=element_blank()) +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5))
dev.off()

## look for overlap
#remove missing
new <- lapply(d, function(x) {
  x <- x[which(x$pi != "na"),]
  x      
})

####### merge all sets to include only overlapping ########

one <- (merge(new[['AAAA_F1_REP1']], new[['AAAA_F1_REP2']], by="id"))[,-c(7,8)]
colnames(one) <- c("id", "gene", "position", "variants_AAAA_F1_REP1", "prop_covered_AAAA_F1_REP1", "pi_AAAA_F1_REP1",
                    "variants_AAAA_F1_REP2", "prop_covered_AAAA_F1_REP2", "pi_AAAA_F1_REP2")

two <- (merge(one, new[['AAAA_F1_REP3']], by="id"))[,-c(10,11)]
colnames(two) <- c(colnames(one), "variants_AAAA_F1_REP3", "prop_covered_AAAA_F1_REP3", "pi_AAAA_F1_REP3")
three <- (merge(two, new[['AAAA_F1_REP4']], by="id"))[,-c(13,14)]
colnames(three) <- c(colnames(two), "variants_AAAA_F1_REP4", "prop_covered_AAAA_F1_REP4", "pi_AAAA_F1_REP4")

four <- (merge(three, new[['HHHH_F1_REP1']], by="id"))[,-c(16,17)]
colnames(four) <- c(colnames(three), "variants_HHHH_F1_REP1", "prop_covered_HHHH_F1_REP1", "pi_HHHH_F1_REP1")

five <- (merge(four, new[['HHHH_F1_REP2']], by="id"))[,-c(19,20)]
colnames(five) <- c(colnames(four), "variants_HHHH_F1_REP2", "prop_covered_HHHH_F1_REP2", "pi_HHHH_F1_REP2")

six <- (merge(five, new[['HHHH_F1_REP3']], by="id"))[,-c(22,23)]
colnames(six) <- c(colnames(five), "variants_HHHH_F1_REP3", "prop_covered_HHHH_F1_REP3", "pi_HHHH_F1_REP3")

df.all <- (merge(six, new[['HHHH_F1_REP4']], by="id"))[,-c(25,26)]
colnames(df.all) <- c(colnames(six), "variants_HHHH_F1_REP4", "prop_covered_HHHH_F1_REP4", "pi_HHHH_F1_REP4")


# calc cutoff values

AAAA_REP1 <- which(df.all$pi_AAAA_F1_REP1 < quantile(df.all$pi_AAAA_F1_REP1, na.rm=TRUE, 0.01))
AAAA_REP2 <- which(df.all$pi_AAAA_F1_REP2 < quantile(df.all$pi_AAAA_F1_REP2, na.rm=TRUE, 0.01))
AAAA_REP3 <- which(df.all$pi_AAAA_F1_REP3 < quantile(df.all$pi_AAAA_F1_REP3, na.rm=TRUE, 0.01))
AAAA_REP4 <- which(df.all$pi_AAAA_F1_REP4 < quantile(df.all$pi_AAAA_F1_REP4, na.rm=TRUE, 0.01))
HHHH_REP1 <- which(df.all$pi_HHHH_F1_REP1 < quantile(df.all$pi_HHHH_F1_REP1, na.rm=TRUE, 0.01))
HHHH_REP2 <- which(df.all$pi_HHHH_F1_REP2 < quantile(df.all$pi_HHHH_F1_REP2, na.rm=TRUE, 0.01))
HHHH_REP3 <- which(df.all$pi_HHHH_F1_REP3 < quantile(df.all$pi_HHHH_F1_REP3, na.rm=TRUE, 0.01))
HHHH_REP4 <- which(df.all$pi_HHHH_F1_REP4 < quantile(df.all$pi_HHHH_F1_REP4, na.rm=TRUE, 0.01))


##
## look for number of low het genes overlapping between 2, 3, and 4 for each group
## 

##### 2 reps
two.1 <- intersect(AAAA_REP1, AAAA_REP2)
two.2 <- intersect(AAAA_REP1, AAAA_REP3)
two.3 <- intersect(AAAA_REP1, AAAA_REP4)
two.4 <- intersect(AAAA_REP2, AAAA_REP3)
two.5 <- intersect(AAAA_REP2, AAAA_REP4)
two.6 <- intersect(AAAA_REP3, AAAA_REP4)
length(unique(c(two.1, two.2, two.3, two.4, two.5, two.6)))
# 63

two.1 <- intersect(HHHH_REP1, HHHH_REP2)
two.2 <- intersect(HHHH_REP1, HHHH_REP3)
two.3 <- intersect(HHHH_REP1, HHHH_REP4)
two.4 <- intersect(HHHH_REP2, HHHH_REP3)
two.5 <- intersect(HHHH_REP2, HHHH_REP4)
two.6 <- intersect(HHHH_REP3, HHHH_REP4)
length(unique(c(two.1, two.2, two.3, two.4, two.5, two.6)))
# 58

##### 3 reps
three.1 <- Reduce(intersect, list(AAAA_REP1, AAAA_REP2, AAAA_REP3))
three.2 <- Reduce(intersect, list(AAAA_REP1, AAAA_REP2, AAAA_REP4))
three.3 <- Reduce(intersect, list(AAAA_REP1, AAAA_REP3, AAAA_REP4))
three.4 <- Reduce(intersect, list(AAAA_REP4, AAAA_REP2, AAAA_REP3))
length(unique(c(three.1, three.2, three.3, three.4)))
# 36

three.1 <- Reduce(intersect, list(HHHH_REP1, HHHH_REP2, HHHH_REP3))
three.2 <- Reduce(intersect, list(HHHH_REP1, HHHH_REP2, HHHH_REP4))
three.3 <- Reduce(intersect, list(HHHH_REP1, HHHH_REP3, HHHH_REP4))
three.4 <- Reduce(intersect, list(HHHH_REP4, HHHH_REP2, HHHH_REP3))
length(unique(c(three.1, three.2, three.3, three.4)))
# 40

##### 4 reps

AAAA_overlap <- Reduce(intersect, list(AAAA_REP1, AAAA_REP2, AAAA_REP3, AAAA_REP4))
HHHH_overlap <- Reduce(intersect, list(HHHH_REP1, HHHH_REP2, HHHH_REP3, HHHH_REP4))
length(AAAA_overlap)
# 21
length(HHHH_overlap)
# 30


# all
length(intersect(Reduce(intersect, list(HHHH_REP1, HHHH_REP2, HHHH_REP3, HHHH_REP4)), 
    Reduce(intersect, list(AAAA_REP1, AAAA_REP2, AAAA_REP3, AAAA_REP4))))

