
## Plotting interesting candidate individual genes.


First plot NADH histogram figure. Fig. S3.  

Then indiv genes, fig S2

```r

library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(dplyr)
library(Rmisc)
library(reshape2)
library(ggthemes)
library(ggpubr)
library(matrixStats)

#dat <- read.csv("~/reciprocal_t/analysis/pi.overlp.bed", header=FALSE, sep="\t")

#colnames(dat) <- c("CHR", "start", "stop", "cov", "nsnp", "pi",
#                    "pi_snp", "gp", "gene", "snp_start", "snp_stop", "sig")

all_out <- read.table("~/reciprocal_t/analysis/AF_change.txt", header=TRUE)


# make histogram of af changes

all_out$f1 <- abs(all_out$HHHH_F1_mean - all_out$AAAA_F1_mean)
all_out$f2 <- abs(all_out$HHHH_F2_mean - all_out$AAAA_F2_mean)
all_out$f3 <- abs(all_out$HHHH_F3_mean - all_out$AAAA_F3_mean)

#png(file="~/reciprocal_t/figures/af_histogram.png", height=7, w=5, units="in", res=300)

#par(mfrow=c(3,1),
#mar = c(4, 4.3, 1, 0)) # Set the margin on all sides to 2

p1 <- ggplot(all_out, aes(x=as.numeric(f1))) + 
  geom_histogram(breaks=seq(0,max(all_out$f1)+ 0.01,0.02),
  color="black", fill="gray55") +
  theme_classic() +
  geom_rug(length = unit(0.05, "npc"), alpha = 0.2) +
   scale_y_continuous(expand = c(0.1, 0.1),
        breaks=c(0,25000, 50000,75000,100000, 125000),
        labels=c("0", "25", "50", "75", "100", "125"),
        limits=c(0,126000),
        name="count (thousands)") +
   scale_x_continuous(expand = c(0.01, 0), limits=c(0, 0.9)) +
   ggtitle("F1") +
   xlab("") +
     theme(plot.title = element_text(hjust = 0.5)) +
   geom_segment(x = mean(all_out$f1), xend = mean(all_out$f1),
                y = 0, yend= 120000, 
                linetype="dashed", color="dodgerblue3", size=1.5) +
   geom_text(x=mean(all_out$f1), label="mean", y=128000, color="black", size= 4)+
   geom_segment(x= max(all_out$f1), xend = max(all_out$f1),
                y = 120000, yend= 0,
                arrow = arrow(length = unit(0.2, "cm"))) +
   geom_text(x=max(all_out$f1), label="NADH", y=128000, color="black", size= 4)

p2 <- ggplot(all_out, aes(x=as.numeric(f2))) + 
  geom_histogram(breaks=seq(0,max(all_out$f2)+ 0.01,0.02),
  color="black", fill="gray55") +
  theme_classic() +
  geom_rug(length = unit(0.05, "npc"), alpha = 0.2) +
   scale_y_continuous(expand = c(0.1, 0.1),
        breaks=c(0,25000, 50000,75000,100000, 125000),
        labels=c("0", "25", "50", "75", "100", "125"),
        name="count (thousands)", 
        limits=c(0,126000)) +
   scale_x_continuous(expand = c(0.01, 0), limits=c(0, 0.9)) +
   ggtitle("F2") +
   xlab("") +
     theme(plot.title = element_text(hjust = 0.5)) +
   geom_segment(x = mean(all_out$f2), xend = mean(all_out$f2),
                y = 0, yend= 120000, 
                linetype="dashed", color="dodgerblue3", size=1.5) +
   geom_text(x=mean(all_out$f2), label="mean", y=128000, color="black", size= 4)+
   geom_segment(x= max(all_out$f2), xend = max(all_out$f2),
                y = 120000, yend= 0,
                arrow = arrow(length = unit(0.2, "cm"))) +
   geom_text(x=max(all_out$f2), label="NADH", y=128000, color="black", size= 4)

p3 <- ggplot(all_out, aes(x=as.numeric(f3))) + 
  geom_histogram(breaks=seq(0,max(all_out$f3)+ 0.01,0.02),
  color="black", fill="gray55") +
  theme_classic() +
  geom_rug(length = unit(0.05, "npc"), alpha = 0.2) +
   scale_y_continuous(expand = c(0.1, 0.1),
        breaks=c(0,25000, 50000,75000,100000, 125000),
        labels=c("0", "25", "50", "75", "100", "125"),
        name="count (thousands)", 
        limits=c(0,126000)) +
   scale_x_continuous(expand = c(0.01, 0), limits=c(0, 0.9)) +
   ggtitle("F3") +
   xlab("Abs difference in allele frequency: AM - OWA") +
     theme(plot.title = element_text(hjust = 0.5)) +
   geom_segment(x = mean(all_out$f3), xend = mean(all_out$f3),
                y = 0, yend= 120000, 
                linetype="dashed", color="dodgerblue3", size=1.5) +
   geom_text(x=mean(all_out$f3), label="mean", y=128000, color="black", size= 4)+
   geom_segment(x= max(all_out$f3), xend = max(all_out$f3),
                y = 120000, yend= 0,
                arrow = arrow(length = unit(0.2, "cm"))) +
   geom_text(x=max(all_out$f3), label="NADH", y=128000, color="black", size= 4)

ggsave(filename="~/reciprocal_t/figures/af_histogram.png", ggarrange(p1, p2, p3, nrow=3, ncol=1),
    h=7, w=5, units="in", dpi=300)



mean(c(all_out$f1[which(all_out$f1 > 0.65)],all_out$f2[which(all_out$f2 > 0.7)],all_out$f3[which(all_out$f3 > 0.75)]))


####################################################################################
##########################################
# make plots of indiv genes
##########################################
####################################################################################


all <- all_out

# find SE of delta af:

all$rep1_f1 <- all$AAAA_F1_REP1 - all$HHHH_F1_REP1
all$rep2_f1 <- all$AAAA_F1_REP2 - all$HHHH_F1_REP2
all$rep3_f1 <- all$AAAA_F1_REP3 - all$HHHH_F1_REP3
all$rep4_f1 <- all$AAAA_F1_REP4 - all$HHHH_F1_REP4

all$rep1_f2 <- all$AAAA_F2_REP1 - all$HHHH_F2_REP1
all$rep2_f2 <- all$AAAA_F2_REP2 - all$HHHH_F2_REP2
all$rep3_f2 <- all$AAAA_F2_REP3 - all$HHHH_F2_REP3
all$rep4_f2 <- all$AAAA_F2_REP4 - all$HHHH_F2_REP4

all$rep1_f3 <- all$AAAA_F3_REP1 - all$HHHH_F3_REP1
all$rep2_f3 <- all$AAAA_F3_REP2 - all$HHHH_F3_REP2
all$rep3_f3 <- all$AAAA_F3_REP3 - all$HHHH_F3_REP3
all$rep4_f3 <- all$AAAA_F3_REP4 - all$HHHH_F3_REP4

all$rep1_af <- rowMeans((all[,grep("rep1_f1|rep1_f2|rep1_f3", colnames(all))]))
all$rep2_af <- rowMeans((all[,grep("rep2_f1|rep2_f2|rep2_f3", colnames(all))]))
all$rep3_af <- rowMeans((all[,grep("rep3_f1|rep3_f2|rep3_f3", colnames(all))]))
all$rep4_af <- rowMeans((all[,grep("rep4_f1|rep4_f2|rep4_f3", colnames(all))]))


all$af_rep_mean <- rowMeans((all[,grep("rep1_af|rep2_af|rep3_af|rep4_af", colnames(all))]))
# find sd:
all$af_rep_sd <- rowSds(as.matrix(all[,grep("rep1_af|rep2_af|rep3_af|rep4_af", colnames(all))]))
all$af_rep_se <- all$af_rep_sd/sqrt(4)


#####################################
# find genes with large changes
#####################################


large <- all[which(abs(all$af_rep_mean) > quantile(abs(all$af_rep_mean), 0.90)),]
large$gene <- (str_split_fixed(large$SNP, ":", n=2))[,1]

length(unique(large$gene))

# pull out the gene names

write.table(file="~/reciprocal_t/analysis/orf/large_change.txt", unique(large$gene), 
            col.names=F, row.names=FALSE, quote=FALSE, sep="\t")



```


### find ORFs

```bash

cat /data/copepods/tonsa_transcriptome/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.fasta | grep -A1 'TRINITY_DN146090_c0_g4' > ~/reciprocal_t/analysis/orf/tonsa_for_orf.fasta
cat /data/copepods/tonsa_transcriptome/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.fasta | grep -A1 'TRINITY_DN150080_c2_g2' >> ~/reciprocal_t/analysis/orf/tonsa_for_orf.fasta
cat /data/copepods/tonsa_transcriptome/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.fasta | grep -A1 'TRINITY_DN148665_c4_g2' >> ~/reciprocal_t/analysis/orf/tonsa_for_orf.fasta
cat /data/copepods/tonsa_transcriptome/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.fasta | grep -A1 'TRINITY_DN141476_c0_g3' >> ~/reciprocal_t/analysis/orf/tonsa_for_orf.fasta

cat /data/copepods/tonsa_transcriptome/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.fasta |\
 grep -A1 -f ~/reciprocal_t/analysis/orf/large_change.txt >> ~/reciprocal_t/analysis/orf/tonsa_for_orf.fasta

sed '/--/d' tonsa_for_orf.fasta > tonsa_for_orf.fasta1
mv tonsa_for_orf.fasta1 tonsa_for_orf.fasta

~/Documents/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -t ~/Documents/UVM/Reciprocal_transplant/orf/tonsa_for_orf.fasta

~/Documents/TransDecoder-TransDecoder-v5.5.0/TransDecoder.Predict -t ~/Documents/UVM/Reciprocal_transplant/orf/tonsa_for_orf.fasta


# move fasta file and gff to snp location:

cp ~/reciprocal_t/analysis/orf/tonsa_for_orf.fasta ~/bin/snpEff/data/tonsa/sequences.fa

#scp ~/Documents/UVM/Reciprocal_transplant/orf/tonsa_for_orf.fasta.transdecoder.gff3  'rbrennan@pespenilab.uvm.edu:~/bin/snpEff/data/tonsa'

mv ~/bin/snpEff/data/tonsa/tonsa_for_orf.fasta.transdecoder.gff3 genes.gff

cd ~/bin/snpEff/data/tonsa

# build the snp eff files
java -jar ~/bin/snpEff/snpEff.jar build -gff3 -v tonsa

#make vcf for the genes of interest:

#cat filtered_variants.sync| grep -F 'TRINITY_DN146090_c0_g4|TRINITY_DN150080_c2_g2|TRINITY_DN148665_c4_g2|TRINITY_DN141476_c0_g3'  > ~/reciprocal_t/analysis/orf/selected_variants.txt

cat ~/reciprocal_t/analysis/filtered_variants.sync | grep -f ~/reciprocal_t/analysis/orf/large_change.txt > ~/reciprocal_t/analysis/orf/selected_variants.txt


```


```r

# make vcf for snp eff to use
library(poolfstat)

snps <- popsync2pooldata(sync.file = "~/reciprocal_t/analysis/orf/selected_variants.txt", 
    poolsizes=rep(50,48), min.maf=0.0001)

df <- as.data.frame(snps@snp.info)
df$SNP <- paste(df$V1, df$V2, sep=":")
df <- df[df$SNP %in% large$SNP,]

out <- data.frame(CHROM = df$V1, POS=df$V2, ID=paste(df$V1, df$V2, sep=":"),
                REF=df$V3, ALT = df$V4, QUAL = rep(40, nrow(df)), FILTER=rep("PASS", nrow(df)),
                INFO = rep("DP=50",nrow(df)), FORMAT = rep("GT:DP",nrow(df)),
                SAMP1 = rep("0/1:50",nrow(df)))


write.table(file="~/reciprocal_t/analysis/orf/selected_variants.vcf", out, 
            col.names=T, row.names=FALSE, quote=FALSE, sep="\t")

```


run snpeff to get effect of snps

```bash
java -Xmx4g -jar ~/bin/snpEff/snpEff.jar -c ~/bin/snpEff/snpEff.config  \
-v tonsa  ~/reciprocal_t/analysis/orf/selected_variants.vcf  \
> ~/reciprocal_t/analysis/orf/selected_variants_ann.vcf > ~/reciprocal_t/analysis/orf/snpEff.stdout 2> ~/reciprocal_t/analysis/orf/snpEff.stderr
```


process the variant effects

```r

library(stringr)
library(ggplot2)
library(gridExtra)
library(MASS)

dat <- read.table("~/reciprocal_t/analysis/orf/selected_variants_ann.vcf", stringsAsFactors=FALSE)

# split annotations
eff <- strsplit(as.character(dat$V8), split=",", fixed=TRUE)

out <- dat

for (i in 1:length(eff)){
    if(length(eff[[i]]) == 1){
        out$V8[i] <- eff[[i]]
    }
    else{
        # pull out the effect of each variant
        # they're ordered by severity
        out$V8[i] <- eff[[i]][1]
        }
    }

# add snp id
out$SNP <- paste(out$V1, out$V2, sep=":")

out_ann <- data.frame(SNP=out$SNP, ANN=out$V8)
# split annotation
ann_split <- (str_split_fixed(out$V8, "\\|", n=16))
out_ann <-(cbind(out$SNP, ann_split))

colnames(out_ann) <- c("SNP","Allele","Annotation","Annotation_Impact","Gene_Name","Gene_ID", "Feature_Type","Feature_ID", "Transcript_BioType","Rank","HGVS.c","HGVS.p","cDNA.pos",  "cDNA.length","CDS.pos","Distance","WARNINGS")

########################################
# assign functional categories
########################################

new <- as.data.frame(out_ann)

# assign class  ie, intergenic, etc
new$class <- c(NA)

new$class[which(new$Annotation == "downstream_gene_variant" |
    new$Annotation == "intergenic_region")] <- c("intergenic")
new$class[which(new$Annotation == "upstream_gene_variant")] <- c("upstream")
new$class[which(new$Annotation == "intron_variant" |
    new$Annotation == "splice_region_variant&intron_variant"|
    new$Annotation == "splice_acceptor_variant&intron_variant"|
    new$Annotation == "splice_donor_variant&splice_region_variant&intron_variant"|
    new$Annotation == "splice_region_variant"|
    new$Annotation == "splice_region_variant&intron_variant" |
    new$Annotation == "splice_region_variant&synonymous_variant"|
    new$Annotation == "splice_region_variant&non_coding_transcript_exon_variant"|
    new$Annotation == "non_coding_transcript_exon_variant"|
    new$Annotation == "splice_donor_variant&intron_variant"|
    new$Annotation == "3_prime_UTR_variant"|
    new$Annotation == "splice_acceptor_variant&splice_donor_variant&intron_variant"|
    new$Annotation == "splice_acceptor_variant&splice_region_variant&intron_variant")] <- c("intron")
new$class[which(new$Annotation == "synonymous_variant" |
    new$Annotation == "stop_retained_variant"|
    new$Annotation == "splice_region_variant&stop_retained_variant")] <- c("synonymous")
new$class[which(new$Annotation == "missense_variant" |
    new$Annotation == "stop_gained" |
    new$Annotation == "stop_lost"|
    new$Annotation == "missense_variant&splice_region_variant"|
    new$Annotation == "start_lost"|
    new$Annotation == "stop_gained&splice_region_variant"|
    new$Annotation == "stop_lost&splice_region_variant"|
    new$Annotation == "frameshift_variant"|
    new$Annotation == "disruptive_inframe_deletion"|
    new$Annotation == "frameshift_variant&stop_gained"|
    new$Annotation == "conservative_inframe_insertion"|
    new$Annotation == "conservative_inframe_deletion"|
    new$Annotation == "start_lost&conservative_inframe_deletion"|
    new$Annotation == "disruptive_inframe_insertion"|
    new$Annotation == "stop_gained&disruptive_inframe_insertion"|
    new$Annotation == "frameshift_variant&splice_region_variant"|
    new$Annotation == "frameshift_variant&stop_lost"|
    new$Annotation == "frameshift_variant&splice_region_variant"|
    new$Annotation == "frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant"|
    new$Annotation == "initiator_codon_variant")] <- c("non-synonymous")


new$gene <- (str_split_fixed(out$SNP, ":", n=2))[,1]
large$gene <- (str_split_fixed(large$SNP, ":", n=2))[,1]
new.sub <- new[new$SNP %in% large$SNP,]
# now 178
new <- new.sub

ns <- new[which(new$class == "non-synonymous"),]

ns$gene <- str_split_fixed(ns$SNP, ":", n=2)[,1]

write.table(file="~/reciprocal_t/analysis/orf/ns_change.txt", unique(ns$gene), 
            col.names=F, row.names=FALSE, quote=FALSE, sep="\t")


```


I manually went through the ns_change.txt file above to find candidates. for example w/ code below.


```bash
zcat /data/copepods/tonsa_transcriptome/Atonsa_gen_trans_agp_gff/blast2go_gff_export_20170613_0815.gff.gz |\
 grep -f ~/reciprocal_t/analysis/orf/ns_change.txt | grep 'ATPase'

```


plot the interesting genes.

```r

all$gene <- (str_split_fixed(all$SNP, ":", n=2))[,1]

gene_in <- data.frame(
    gene = c("TRINITY_DN141476_c0_g3","TRINITY_DN150080_c2_g2","TRINITY_DN148665_c4_g2",
        "TRINITY_DN147294_c5_g2", "TRINITY_DN149301_c0_g7", "TRINITY_DN149639_c2_g19",
             "TRINITY_DN147803_c0_g1", "TRINITY_DN151143_c1_g7", 
             "TRINITY_DN151056_c4_g4", "TRINITY_DN119996_c0_g1", "TRINITY_DN141181_c1_g5",
             "TRINITY_DN144102_c0_g16", "TRINITY_DN151024_c6_g2", "TRINITY_DN149855_c0_g15",
             "TRINITY_DN146090_c0_g4"),
    id = c("NADH-ubiquinone oxidoreductase 49 kDa subunit",
        "vitellogenin receptor",
        "ribosome biogenesis BMS1",
        "activator of 90 kDa heat shock ATPase",
        "heat shock 75 mitochondrial",
        "methyl CpG binding protein 2",
        "actin cytoskeleton-regulatory complex PAN1",
        "Transcription initiation factor TFIID subunit 2",
        "cation-transporting ATPase 13A3",
        "ubiquinol-cytochrome-c reductase complex assembly factor 2",
        "general transcription factor IIF subunit 1",
        "ribosome biogenesis BOP1 homolog",
        "round spermatid basic 1",
        "vitellogenin 2",
        "ribosomal S6 kinase"
        )
    )

mdat <- merge(all, gene_in, by="gene")

mdat$POS <- as.numeric(str_split_fixed(mdat$SNP, ":", n=2)[,2])

# which points are non-syn?
# mdat$Non_syn <- c("synonymous")

sub1 <- subset(mdat, gene == "TRINITY_DN119996_c0_g1")
sub2 <- subset(mdat, gene == "TRINITY_DN141181_c1_g5")
sub3 <- subset(mdat, gene == "TRINITY_DN141476_c0_g3" & POS < 550)
sub4 <- subset(mdat, gene == "TRINITY_DN144102_c0_g16" & POS < 1500 & POS > 500)
sub5 <- subset(mdat, gene == "TRINITY_DN147294_c5_g2")
sub6 <- subset(mdat, gene == "TRINITY_DN147803_c0_g1" & POS < 5000 & POS > 2000)
sub7 <- subset(mdat, gene == "TRINITY_DN148665_c4_g2" & POS < 2000 & POS > 1000)
sub8 <- subset(mdat, gene == "TRINITY_DN149301_c0_g7")
sub9 <- subset(mdat, gene == "TRINITY_DN149639_c2_g19")
sub10 <- subset(mdat, gene == "TRINITY_DN149855_c0_g15")
sub11 <- subset(mdat, gene == "TRINITY_DN150080_c2_g2" & POS < 6000 & POS > 5000)
sub12 <- subset(mdat, gene == "TRINITY_DN151024_c6_g2" & POS < 3800)
sub13 <- subset(mdat, gene == "TRINITY_DN151056_c4_g4" & POS < 2500)
sub14 <- subset(mdat, gene == "TRINITY_DN151143_c1_g7")
sub15 <- subset(mdat, gene == "TRINITY_DN146090_c0_g4")


submdat <- rbind(sub1, sub2, sub3, sub4, sub5, sub6, sub7, sub8, sub9, sub10, 
                sub11, sub12, sub13, sub14, sub15)

library(stringr)
var_width = 30
my_plot_data <- mutate(submdat, pretty_varname = str_wrap(id, width = var_width))

sig_pt <- my_plot_data[my_plot_data$SNP %in% ns$SNP,]

p <- ggplot(my_plot_data, aes(x=POS, y=abs(af_rep_mean))) +
#   geom_line(lty=3, lwd=0.8, size=05) +
    geom_errorbar(aes(ymin=abs(af_rep_mean)-af_rep_se, ymax=abs(af_rep_mean)+af_rep_se), 
        colour="black", width=.1) +
    geom_point(shape=21, fill="gray55") +
    facet_wrap( ~ pretty_varname, scales = "free", ncol=3) +
    theme_bw() +
    geom_point(data=sig_pt, aes(x=POS, y=abs(af_rep_mean)), shape=21, fill=c("darkorange1"), size=1.8)+
    ylab("Change in allele frequency\n OWA vs AM") +
    xlab("Position") +
    theme(strip.text = element_text(size=6)) +
    geom_hline(yintercept=mean(abs(all$af_rep_mean)), linetype="solid", 
        color="dodgerblue3", size=0.7)
p

ggsave(file="~/reciprocal_t/figures/af_scatters.png", p, h=8, w=7)


```
