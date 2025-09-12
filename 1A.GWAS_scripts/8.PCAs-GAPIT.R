library("readxl")
library("RColorBrewer")
library("ggplot2")

library(dplyr)
library(tidyr)
library(stringr)
library(viridis)
library(plyr)
library(rcartocolor)
library(ggpubr)
#install.packages("rgl")
library(rgl)

rm(list = ls())


##Colour palettes 
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
safe_pal_5 <- c("#000000", "#CC6677", "#CC5500", "#117733", "#6699CC")
safe_pal_7 <- c("#000000", "#CC6677", "#5EC962", "#117733", "#332288", "#6699CC", "#CC5500" )
safe_pal_8 <- c("#000000", "#888888", "#5EC962", "#117733", "#332288", "#6699CC", "#CC5500","#AA4499")

safe_pal_8 <- c("#117733","#000000", "#888888", "#332288", "#5EC962","#AA4499", "#6699CC","#CC5500")

safe_pal_9 <- c("#117733", "#990000", "#888888", "#000000", "#332288", "#5EC962", "#AA4499", "#6699CC", "#CC5500")

setwd("/STATS/23.02-GWAS/pca3_thin5_bytrait/")

GAPIT.Genotype.PCA <- read.csv("GAPIT.Genotype.PCA.csv", header=TRUE)
colnames(GAPIT.Genotype.PCA) <- c("TAXA", "PC1", "PC2", "PC3")
GAPIT_EV <- read.csv("GAPIT.Genotype.PCA_eigenvalues.csv", header=TRUE)
#Seed_phenotypes<- read_excel("Bean project/Common bean samples/22.02 Seed phenotypes/22.02-updated-phenotypes-R.xlsx")

K6_clustering <- read_excel("Bean project/STATS/22.01 Variant calling - GATK/22.11-ADMIXTURE-LD0.5/Q_files/K6_Q_files/K6_clustering_PCA-phenotypes over 0.7 ADMX.xlsx")

#GAPIT_pca_PHENOS <- merge(GAPIT.Genotype.PCA, Seed_phenotypes,by = "Sample_name")
GAPIT_pca_PHENOS_ADMX <- merge(GAPIT.Genotype.PCA, K6_clustering, by = "TAXA" )

phenos <- read_xlsx("../23.06-phenotypes-GAPIT.xlsx")
Det <- phenos[, c(1, 4)]

GAPIT_pca_PHENOS_ADMX <- merge(GAPIT_pca_PHENOS_ADMX, Det, by = "TAXA" )
phenos <- merge(GAPIT_pca_PHENOS_ADMX, phenos, by = "TAXA" )

### T-test
phenos$pho
test <- t.test(Admixed_K6_C ~ `photo_sens`, data = phenos, alternative = "two.sided", var.equal = FALSE)



###### Plots ############## ------------------------------------------
#Shape countries
#Posters size=4

pdf("PCAs/PCA_byphoto_sensandcluster_K6_AM-match-admx-Col_GAPIT_poster.pdf", height=6, width = 7)
p1_2_PI <- ggplot(GAPIT_pca_PHENOS_ADMX,  aes(x= PC1 , y= PC2)) +
  geom_point() + theme_classic() + ylab("PC 2 (5.06%)") + #PC 2 (5.06%) PC 3 (3.67%)
  xlab("PCA 1 (38.8%)") + geom_point(size = 2.5, stroke = 0.5) + 
  aes(colour = `photo sens`, stroke = 1) + #colour = Admixed_K6_C_K2
  labs(colour = "Photoperiod \nSensitivity") + #, colour="Ancestry"
  scale_color_manual(values = c2, labels = c("Photoperiod \ninsensitive", "Photoperiod \nsensitive")) + 
  theme(text = element_text(size=12), legend.title=element_text(size=12)) 
  #scale_shape_manual(labels = c("Photoperiod \ninsensitive", "Photoperiod \nsensitive"), values = c(15:16))
  #guides(color = "none") + 
#+ scale_shape_manual(values=1:10)
#+ geom_label_repel(aes(label=pca_0.5_merge_K6$'Tag 2'), max.overlaps = Inf)
p1_2_PI
dev.off()
safe_pal_8 <- c("#117733", "#000000", "#888888", "#332288", "#5EC962", "#AA4499", "#6699CC", "#CC5500")

c2 <- c(
  "darkorange4", # purple
 "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
p1_2_SS <- ggplot(GAPIT_pca_PHENOS_ADMX,  aes(x= PC1 , y= PC2)) +
  geom_point() + theme_classic() + ylab("PC 2 (5.06%)") + #PC 2 (5.06%) PC 3 (3.67%)
  xlab("PCA 1 (38.8%)") + geom_point(size = 2.5, stroke = 0.5) + 
  aes(colour = `Estimated seed size`, stroke = 1) + #colour = Admixed_K6_C_K2,
  labs(colour = "Estimated \nSeed Size") + #, colour="Ancestry"
  scale_color_manual(values = c2) + 
  theme(text = element_text(size=12), legend.title=element_text(size=12)) 
  #guides(color = "none")
#+ scale_shape_manual(values=1:10)
#+ geom_label_repel(aes(label=pca_0.5_merge_K6$'Tag 2'), max.overlaps = Inf)
p1_2_SS

GAPIT_pca_PHENOS_ADMX$Country[GAPIT_pca_PHENOS_ADMX$Country == "Heirlooms(NA)"] <- "Heirlooms"
c25 <- c(
  "dodgerblue2", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "orchid1", # lt pink
  "lightgreen",
  "maroon", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

p1_2_C <- ggplot(GAPIT_pca_PHENOS_ADMX,  aes(x= PC1 , y= PC2)) +
  geom_point() + theme_classic() + ylab("PC 2 (5.06%)") + #PC 2 (5.06%) PC 3 (3.67%)
  xlab("PCA 1 (38.8%)") + geom_point(size = 2.5, stroke = 0.5) + 
  aes(colour = Country, stroke = 1) + #, colour = Admixed_K6_C_K2
  labs(colour="Country") + scale_color_manual(values = c25) + #, colour="Ancestry"
  theme(text = element_text(size=12), legend.title=element_text(size=12)) + 
  scale_shape_manual(values=c(0:4, 8, 11, 12, 15:17)) 
  #guides(color = "none")
#+ geom_label_repel(aes(label=pca_0.5_merge_K6$'Tag 2'), max.overlaps = Inf)
p1_2_C

GAPIT_pca_PHENOS_ADMX$Determinancy
GAPIT_pca_PHENOS_ADMX$Determinancy[GAPIT_pca_PHENOS_ADMX$Determinancy == "0"] <- "Determinate"
GAPIT_pca_PHENOS_ADMX$Determinancy[GAPIT_pca_PHENOS_ADMX$Determinancy == "1"] <- "Indeterminate"

c25 <- c( 
  "maroon", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
p1_2_D <- ggplot(GAPIT_pca_PHENOS_ADMX,  aes(x= PC1 , y= PC2)) +
  geom_point() + theme_classic() + ylab("PC 2 (5.06%)") + #PC 2 (5.06%) PC 3 (3.67%)
  xlab("PCA 1 (38.8%)") + geom_point(size = 2.5, stroke = 0.5) + 
  aes(colour = Determinancy,  stroke = 1) + #colour = Admixed_K6_C_K2
  labs(colour = "Determinate") + #, colour="Ancestry"
  scale_color_manual(values = c25) + 
  theme(text = element_text(size=12), legend.title=element_text(size=12)) 
  #guides(color = "none")
#+ scale_shape_manual(values=1:10)
#+ geom_label_repel(aes(label=pca_0.5_merge_K6$'Tag 2'), max.overlaps = Inf)
p1_2_D


#K6 
pdf("PCA_cluster_K6_AM-match-admx-Col_GAPIT_1_3.pdf", height=6, width = 7)
p1_3 <- ggplot(GAPIT_pca_PHENOS_ADMX,  aes(x= PC1 , y= PC3)) +
  geom_point() + theme_classic() + ylab("PC 3 (3.67%)") + #PC 2 (5.06%) PC 3 (3.67%)
  xlab("PC 1 (38.8%)") + geom_point(size=2.5, stroke = 0.5) + aes(colour = Admixed_K6_C_K2) +
  labs(colour="Ancestry") + scale_color_manual(values = safe_pal_8) +
  theme(legend.title=element_text(size=12))
  #labs(tag = "A") 
#+ geom_label_repel(aes(label=pca_0.5_merge_K6$'Tag 2'), max.overlaps = Inf)
p1_3
dev.off()

#K6 
pdf("PCA_cluster_K6_AM-match-admx-Col_GAPIT_1_2.pdf", height=6, width = 7)
p1_2 <- ggplot(GAPIT_pca_PHENOS_ADMX,  aes(x= PC1 , y= PC2)) +
  geom_point() + theme_classic() + ylab("PC 2 (5.06%)") + #PC 2 (5.06%) PC 3 (3.67%)
  xlab("PC 1 (38.8%)") + geom_point(size=2.5, stroke = 0.5) + aes(colour = Admixed_K6_C_K2) +
  labs(colour="Ancestry") + scale_color_manual(values = safe_pal_8) +
  theme(legend.title=element_text(size=12))
#labs(tag = "A") 
#+ geom_label_repel(aes(label=pca_0.5_merge_K6$'Tag 2'), max.overlaps = Inf)
p1_2
dev.off()

#Merge with Heterozygosity plot

p5 = ggarrange(p1, A, ncol = 2, nrow = 1, legend = "right", common.legend = T)
p5

pdf(file = "PCA_GAPIT_K6_A_hetero_accession.pdf",width = 10, height = 5)
plot(p5)
dev.off()

blank_plot <- ggplot() + 
  theme_void() + 
  theme(panel.background = element_blank())

p6 <- ggarrange(p1_2, blank_plot, p1_3,
                ncol = 3, nrow = 1, legend = "right", common.legend = T, 
                labels = c("A","", "B"), widths = c(1, 0.3, 1))


p7 <- ggarrange(p1_2_PI, p1_2_D, p1_2_SS, p1_2_C,
                ncol = 2, nrow = 2, legend = "right", 
                labels = c("C", "D", "E", "F"))

p8 <- ggarrange(p6, p7, nrow = 2, ncol = 1, heights = c(1, 2))


pdf(file = "PCA1_3_2_GAPIT_K6_A_hetero_accession_test.pdf",width = 10, height = 5)
plot(p6)
dev.off()

pdf(file = "PCA1_3_2_GAPIT_K6_A_hetero_accession_phenos_col2.pdf",width = 10, height = 9)
plot(p8)
dev.off()

#Andean PCA #############

setwd("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Bean project/STATS/23.02-GWAS/pca3_thin5_Andean/")

GAPIT.Genotype.PCA.A <- read.csv("GAPIT.Genotype.PCA.csv", header=TRUE)
colnames(GAPIT.Genotype.PCA.A) <- c("TAXA", "PC1", "PC2", "PC3")
GAPIT_EV_A <- read.csv("GAPIT.Genotype.PCA_eigenvalues.csv", header=TRUE)

GAPIT_pca_PHENOS_ADMX_A <- merge(GAPIT.Genotype.PCA.A, K6_clustering, by = "TAXA" )

pdf("PCA_cluster_K6_AM-match-admx-Col_GAPIT.pdf", height=6, width = 7)
p2 <- ggplot(GAPIT_pca_PHENOS_ADMX_A,  aes(x= PC1 , y= PC2)) +
  geom_point() + theme_classic() + ylab("PC 2 (6.1%)") +
  xlab("PC 1 (10.3%)") + geom_point(size=2.5) + aes(colour = Admixed_K6_C_K2) +
  labs(colour="Ancestry") + scale_color_manual(values = safe_pal_8)
#+ geom_label_repel(aes(label=pca_0.5_merge_K6$'Tag 2'), max.overlaps = Inf)
p2
dev.off()

#Merge with Heterozygosity plot

p5 = ggarrange(p2, A, ncol = 2, nrow = 1, legend = "right", common.legend = T,
               labels = "AUTO")
p5

p6 = ggarrange(p1, p5, ncol = 2, nrow = 1, legend = "right", common.legend = T)
p6

pdf(file = "PCA_GAPIT_K6_Andean_GAPIT_hetero_accession_test.pdf",width = 9, height = 4)
plot(p5)
dev.off()


################# Plot line graph eigenvectors ------------------------------------------
x <- sum(GAPIT_EV$x)

GAPIT_EV$percent <- (GAPIT_EV$x/x)*100
GAPIT_EV$y <- c(1:144)

ROWS <- head(GAPIT_EV, 10)

pdf("EIGENVECTORS_GAPIT.pdf", height=6, width = 6)
EV_p <- ggplot(data=ROWS, aes(x=y, y=percent)) +
  geom_line() + ylab("Percentage") + xlab("Principal Component") + 
  geom_point() + theme_bw()
EV_p
dev.off()


################# Heatmap and tree ------------------------------------------

GAPIT.Geno.Kin_Zhang <- read.csv("GAPIT.Genotype.Kin_Zhang.csv", header=F)
names(GAPIT.Geno.Kin_Zhang)[names(GAPIT.Geno.Kin_Zhang) == 'V1'] <- 'TAXA'
K6_clus <- K6_clustering[,c(46,2)]

GAPIT_kin_ADMX <- merge(GAPIT.Geno.Kin_Zhang, K6_clus, by = "TAXA" )
GAPIT_kin_ADMX <- GAPIT_kin_ADMX[, c(146,2:144)]


#install.packages("pheatmap")
library(pheatmap)

