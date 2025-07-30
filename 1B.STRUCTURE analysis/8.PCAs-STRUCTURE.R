library('svglite')
library(SNPRelate)
library(gdsfmt)
library("readxl")
library("RColorBrewer")
library("ggplot2")
library(dplyr)
library(tidyr)
library(stringr)
library(ggfortify)
library(ggrepel)
library(pca3d)
library(reshape)
library("GGally") 
library("data.table")
library(paletteer)
library(viridis)
library(plyr)
library("devtools")
#install.packages("rcartocolor")
library(rcartocolor)

#--------------------------------------
#Colour blind friendly palette over 10
safe_pal <- carto_pal(12, "Safe")
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
safe_pal_5 <- c("#000000", "#CC6677", "#CC5500", "#117733", "#6699CC")
safe_pal_7 <- c("#000000", "#CC6677", "#5EC962", "#117733", "#332288", "#6699CC", "#CC5500" )
safe_pal_8 <- c("#000000", "#CC6677", "#5EC962", "#117733", "#332288", "#6699CC", "#CC5500","#AA4499")
scales::show_col(safe_pal_5)
scales::show_col(safe_pal_7)
#--------------------------------------

rm(list = ls())

Seed_phenotypes<- read_excel("Bean project/Common bean samples/22.02 Seed phenotypes/22.02-updated-phenotypes-R.xlsx")

setwd("Bean project/STATS/22.01 Variant calling - GATK/23.01-STRUCTURE/")

K6_clustering <- 
  read.csv(
    "K6_Clusters over 0.7.csv")
names(K6_clustering)[names(K6_clustering) == 'Sample.name'] <- 'Sample_name'


#On files after filtering for max heterzygosity then ...
################ LD 0.5 ###########################
# variance proportion (%)
#pc.percent_0.5 <- pca_0.5$varprop*100
#head(round(pc.percent_0.5, 2))
#11.29  2.29  2.16  1.65  1.52  1.45

pca_0.5 <- 
  read.csv("Bean project/STATS/22.01 Variant calling - GATK/22.04-filtered_GATK/22.06-Linkage_disequilibrium-in-HPC/PCA_0.5/pca_0.5.txt",
           header=TRUE)

pca_0.5_merge <- merge(pca_0.5, Seed_phenotypes,by = "Sample_name")
pca_0.5_merge_K6 <- merge(pca_0.5_merge, K6_clustering, by = "Sample_name")

write.csv(
  file = "K6_clustering_PCA-phenotypes > 0.7.csv", 
  pca_0.5_merge_K6
)

#------------Plots--------------
setwd("/Bean project/STATS/22.01 Variant calling - GATK/22.04-filtered_GATK")


pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5/PCAs-cluster over 0.7-STRUCTURE/PCA_0.5_cluster-colour_K6.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_merge_K6,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (2.29%)") +
  xlab("PCA 1 (11.29%)") + geom_point(size=3) + aes(colour = Clusters) +
  labs(colour = "Clusters") + scale_color_manual(values = safe_pal_7) 
p
dev.off()

#---------------------------------------------------------
###### Swap shapes and colours ###########

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5/PCAs-cluster over 0.7-STRUCTURE/PCA_0.5_bycountryandcluster_K2_tags.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_merge_K2,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (2.29%)") +
  xlab("PCA 1 (11.29%)") + geom_point(size=4) + aes(colour = Clusters, shape = Country, stroke = 1) +
  labs(shape = "Country of Origin", colour="Clusters") + scale_color_manual(values = safe_pal_5) +
  scale_shape_manual(values=1:10) + geom_label_repel(aes(label=pca_0.5_merge_K2$'Tag 3'), max.overlaps = Inf)
p
dev.off()

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5/PCAs-cluster over 0.7-STRUCTURE/PCA_0.5_bycountryandcluster_K6.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_merge_K6,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (2.29%)") +
  xlab("PCA 1 (11.29%)") + geom_point(size=4) + aes(colour = Clusters, shape = Country, stroke = 1) +
  labs(shape = "Country of Origin", colour="Clusters") + scale_color_manual(values = safe_pal_7) +
  scale_shape_manual(values=1:10) #+ geom_label_repel(aes(label=pca_0.5_merge_K6$'Tag 2'), max.overlaps = Inf)
p
dev.off()

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5/PCAs-cluster over 0.7-STRUCTURE/PCA_0.5_bycountryandcluster_K6_1V3.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_merge_K6,  aes(x= EV1 , y= EV3)) +
  geom_point() + theme_classic() + ylab("PCA 3 (2.16%)") +
  xlab("PCA 1 (11.29%)") + geom_point(size=4) + aes(colour = Clusters, shape = Country, stroke = 1) +
  labs(shape = "Country of Origin", colour="Clusters") + scale_color_manual(values = safe_pal_7) +
  scale_shape_manual(values=1:10) #+ geom_label_repel(aes(label=pca_0.5_merge_K6$'Tag 2'), max.overlaps = Inf)
p
dev.off()

#-------------------------------------------------------------------------------------------------------------------
####### LD 0.5 w.10 - thinned, window size 10  ################
# variance proportion (%)
#pc.percent_0.5_w10 <- pca_0.5_w10$varprop*100
#head(round(pc.percent_0.5_w10, 2))
#21.19  4.15  3.12  2.44  2.03  1.78


pca_0.5_w10 <- 
  read.csv("Bean project/STATS/22.01 Variant calling - GATK/22.04-filtered_GATK/22.06-Linkage_disequilibrium-in-HPC/PCA_0.5_thin_w10/pca_0.5_w10.txt",
           header=TRUE)

pca_0.5_w10_merge <- merge(pca_0.5_w10, Seed_phenotypes,by = "Sample_name")
pca_0.5_w10_merge_K6 <- merge(pca_0.5_w10_merge, K6_clustering, by = "Sample_name")

#---------------------------------------------------------
###### Swap shapes and colours ###########

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5_thin_w10/PCAs-cluster over 0.7-STRUCTURE/PCA_0.5_bycountryandcluster_K2_tags.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_w10_merge_K2,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (4.15%)") +
  xlab("PCA 1 (21.19%)") + geom_point(size=4) + aes(colour = Clusters, shape = Country, stroke = 1) +
  labs(shape = "Country of Origin", colour="Clusters") + scale_color_manual(values = safe_pal_5) +
  scale_shape_manual(values=1:10) + geom_label_repel(aes(label=pca_0.5_w10_merge_K2$'Tag 3'), max.overlaps = Inf)
p
dev.off()

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5_thin_w10/PCAs-cluster over 0.7-STRUCTURE/PCA_0.5_bycountryandcluster_K6.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_w10_merge_K6,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (4.15%)") +
  xlab("PCA 1 (21.19%)") + geom_point(size=4) + aes(colour = Clusters, shape = Country, stroke = 1) +
  labs(shape = "Country of Origin", colour="Clusters") + scale_color_manual(values = safe_pal_7) +
  scale_shape_manual(values=1:10) #+ geom_label_repel(aes(label=pca_0.5_w10_merge_K2$'Tag 2'), max.overlaps = Inf)
p
dev.off()

#---------------------------------------------------------
#Plot 1 v 3 

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5_thin_w10/PCAs-cluster over 0.7-STRUCTURE/PCA_0.5_bycountryandcluster_K2_1V3.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_w10_merge_K2,  aes(x= EV1 , y= EV3)) +
  geom_point() + theme_classic() + ylab("PCA 3 (3.12%)") +
  xlab("PCA 1 (21.19%)") + geom_point(size=4) + aes(colour = Clusters, shape = Country, stroke = 1) +
  labs(shape = "Country of Origin", colour="Clusters") + scale_color_manual(values = safe_pal_5) +
  scale_shape_manual(values=1:10) #+ geom_label_repel(aes(label=pca_0.5_w10_merge_K4$'Tag 2'), max.overlaps = Inf)
p
dev.off()

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5_thin_w10/PCAs-cluster over 0.7-STRUCTURE/PCA_0.5_bycountryandcluster_K6_1V3.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_w10_merge_K6,  aes(x= EV1 , y= EV3)) +
  geom_point() + theme_classic() + ylab("PCA 3 (3.12%)") +
  xlab("PCA 1 (21.19%)") + geom_point(size=4) + aes(colour = Clusters, shape = Country, stroke = 1) +
  labs(shape = "Country of Origin", colour="Clusters") + scale_color_manual(values = safe_pal_7) +
  scale_shape_manual(values=1:10) #+ geom_label_repel(aes(label=pca_0.5_w10_merge_K6$'Tag 2'), max.overlaps = Inf)
p
dev.off()



