#install.packages('svglite')
library('svglite')
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SNPRelate")
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
#install.packages("pca3d")
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
library(ggpubr)
library(readxl)

rm(list = ls())

#--------------------------------------
#Colour blind friendly palette over 10
safe_pal <- carto_pal(12, "Safe")
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
safe_pal_5 <- c("#000000", "#CC6677", "#CC5500", "#117733", "#6699CC")
safe_pal_7 <- c("#000000", "#CC6677", "#5EC962", "#117733", "#332288", "#6699CC", "#CC5500" )
safe_pal_8 <- c("#000000", "#888888", "#5EC962", "#117733", "#332288", "#6699CC", "#CC5500","#AA4499")

#safe_pal_8 <- c("#5EC962", "#117733", "#6699CC", "#AA4499", "#000000", "#888888", "#332288",  "#CC5500")

safe_pal_8 <- c("#117733","#000000", "#888888", "#332288", "#5EC962","#AA4499", "#6699CC","#CC5500")


test <- c("#4a31c5")
scales::show_col(test)
scales::show_col(safe_pal_5)
scales::show_col(safe_pal_8)
scales::show_col(safe_colorblind_palette)
#--------------------------------------

setwd("Bean project/STATS/22.01 Variant calling - GATK/22.04-filtered_GATK")


# Seed_phenotypes<- read_excel("Bean project/Common bean samples/22.02 Seed phenotypes/22.02-updated-phenotypes-R.xlsx")
# Seed_phenotypes_CIAT <- Seed_phenotypes[c(3:92),]

K2_clustering <- 
  read.csv(
    "/STATS/22.01 Variant calling - GATK/22.11-ADMIXTURE-LD0.5/Q_files/K2_Q_files/K2_Clusters over 0.7.csv"
  )


setwd("/Bean project/STATS/22.01 Variant calling - GATK/22.11-ADMIXTURE-LD0.5/Q_files")
K6_clustering <-
  read_excel(
    "Bean project/STATS/22.01 Variant calling - GATK/22.11-ADMIXTURE-LD0.5/Q_files/K6_Q_files/K6_clustering_PCA-phenotypes over 0.7 ADMX.xlsx"
  )

#flowering <- read_csv("/Bean project/22.03 GH Norwich/flowering-between-phenotyping-and-seed-multiplication_ADMX.csv")
#K6_clustering <- merge(K6_clustering, flowering, by = "Sample_name")

#names(K6_clustering)[names(K6_clustering) == 'Sample name'] <- 'Sample_name'
#names(K2_clustering)[names(K2_clustering) == 'Sample.name'] <- 'Sample_name'

#------------------------------------------------------------------------------------------------------------------
#On files after filtering for max heterzygosity then ...
################ LD 0.5 ###########################
# variance proportion (%)
#pc.percent_0.5 <- pca_0.5$varprop*100
#head(round(pc.percent_0.5, 2))
#11.29  2.29  2.16  1.65  1.52  1.45

pca_0.5 <-
  read.csv("Bean project/STATS/22.01 Variant calling - GATK/22.04-filtered_GATK/22.06-Linkage_disequilibrium-in-HPC/PCA_0.5/pca_0.5.txt",
           header=TRUE)

# pca_0.5_merge <- merge(pca_0.5, Seed_phenotypes,by = "Sample_name")
# pca_0.5_merge_K2 <- merge(pca_0.5_merge, K2_clustering, by = "Sample_name")
# pca_0.5_merge_K3 <- merge(pca_0.5_merge, K3_clustering, by = "Sample_name")
# pca_0.5_merge_K4 <- merge(pca_0.5_merge, K4_clustering, by = "Sample_name")
# pca_0.5_merge_K6 <- merge(pca_0.5_merge, K6_clustering, by = "Sample_name")
# pca_0.5_merge_K6 <- K6_clustering

write.csv(
  file = "Bean project/STATS/22.01 Variant calling - GATK/22.11-ADMIXTURE-LD0.5/Q_files/K2_Q_files/K2_clustering_PCA-phenotypes over 0.75.csv", 
  pca_0.5_merge_K2
  )

write.csv(
  file = "Bean project/STATS/22.01 Variant calling - GATK/22.11-ADMIXTURE-LD0.5/Q_files/K6_Q_files/K6_clustering_PCA-phenotypes > 0.7.csv", 
  pca_0.5_merge_K6
)

#K4 <- pca_0.5_merge_K4[,c(1,27)]

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5/PCAs-clusters over 0.7/PCA_0.5_cluster-shape_K6.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_merge_K6,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (2.29%)") +
  xlab("PCA 1 (11.29%)") + geom_point(size=3) + aes(shape = Clusters, stroke = 0.5) +
  labs(shape = "Clusters") + scale_shape_manual(values=1:7)
p
dev.off()

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5/PCAs-clusters over 0.7/PCA_0.5_cluster-colour_K6.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_merge_K6,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (2.29%)") +
  xlab("PCA 1 (11.29%)") + geom_point(size=3) + aes(colour = Clusters) +
  labs(colour = "Clusters") + scale_color_manual(values = safe_pal_7) 
p
dev.off()


pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5/PCAs-clusters > 0.7/PCA_0.5_bycountryandcluster_K3.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_merge_K2,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (2.29%)") +
  xlab("PCA 1 (11.29%)") + geom_point(size=3) + aes(colour = Country, shape = Clusters) +
  labs(colour="Country of Origin", shape = "Clusters") + scale_color_viridis(discrete=TRUE)
p
dev.off()

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5/PCAs-clusters > 0.7/PCA_0.5_bycountryandcluster_K6.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_merge_K6,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (2.29%)") +
  xlab("PCA 1 (11.29%)") + geom_point(size=3) + aes(colour = Country, shape = Clusters) +
  labs(colour="Country of Origin", shape = "Clusters") + scale_color_viridis(discrete=TRUE) +
  scale_shape_manual(values=1:7)
p
dev.off()

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5/PCAs-clusters > 0.7/PCA_0.5_bycountryandcluster_K6_tags.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_merge_K6,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (2.29%)") +
  xlab("PCA 1 (11.29%)") + geom_point(size=3) + aes(colour = Country, shape = Clusters) +
  labs(colour="Country of Origin", shape = "Clusters") + scale_color_viridis(discrete=TRUE) +
  scale_shape_manual(values=1:7) + geom_label_repel(aes(label=pca_0.5_merge_K6$'Tag 2'), max.overlaps = Inf)
p
dev.off()

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5/PCAs-clusters > 0.7/PCA_0.5_bycountryandcluster_K3_tags.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_merge_K3,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (2.29%)") +
  xlab("PCA 1 (11.29%)") + geom_point(size=3) + aes(colour = Country, shape = Clusters) +
  labs(colour="Country of Origin", shape = "Clusters") + scale_color_viridis(discrete=TRUE) + 
  geom_label_repel(aes(label=pca_0.5_merge_K3$'Tag 2'), max.overlaps = Inf) #+
  #scale_shape_manual(values=1:7)
p
dev.off()

#---------------------------------------------------------
###### Swap shapes and colours ###########

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5/PCAs-clusters over 0.7/PCA_0.5_bycountryandcluster_K2_tags.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_merge_K2,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (2.29%)") +
  xlab("PCA 1 (11.29%)") + geom_point(size=4) + aes(colour = Clusters, shape = Country, stroke = 1) +
  labs(shape = "Country of Origin", colour="Clusters") + scale_color_manual(values = safe_pal_5) +
  scale_shape_manual(values=1:10) + geom_label_repel(aes(label=pca_0.5_merge_K2$'Tag 3'), max.overlaps = Inf)
p
dev.off()

#Shape countries
pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5/PCAs-clusters over 0.7/PCA_0.5_bycountryandcluster_K6_AM-match-admx-Col2.pdf", height=6, width = 8)
p1 <- ggplot(K6_clustering,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (2.29%)") +
  xlab("PCA 1 (11.29%)") + geom_point(size=4) + aes(colour = Admixed_K6_C, shape = Country, stroke = 1) +
  labs(shape = "Country", colour="Ancestry") + scale_color_manual(values = safe_pal_9) +
  scale_shape_manual(values=1:10) #+ geom_label_repel(aes(label=pca_0.5_merge_K6$'Tag 2'), max.overlaps = Inf)
p1
dev.off()
safe_pal_9 <- c("#117733", "#000000", "#888888", "#805040", "#332288", "#5EC962", "#AA4499", "#6699CC", "#CC5500")

#scales::show_col(safe_pal_8)

#Shape growth habit
p1 <- ggplot(K6_clustering,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (2.29%)") +
  xlab("PCA 1 (11.29%)") + geom_point(size=3) + aes(colour = Admixed_K6_C_K2, shape = `Growth habit CIAT then 2023 GH`, stroke = 1) +
  labs(shape = "Growth habit", colour="Ancestry") + scale_color_manual(values = safe_pal_8) 
  #+ geom_label_repel(aes(label=pca_0.5_merge_K6$'Tag 2'), max.overlaps = Inf)
p1

#K6 
p1 <- ggplot(K6_clustering,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (2.29%)") +
  xlab("PCA 1 (11.29%)") + geom_point(size=3) + aes(colour = Admixed_K6_C_K2) +
  labs(colour="Ancestry") + scale_color_manual(values = safe_pal_8) 
#+ geom_label_repel(aes(label=pca_0.5_merge_K6$'Tag 2'), max.overlaps = Inf)
p1


#---------------------------------------------------------
#1 v 3 plots

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5/PCAs-clusters over 0.7/PCA_0.5_bycountryandcluster_K2_1V3_tags.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_merge_K2,  aes(x= EV1 , y= EV3)) +
  geom_point() + theme_classic() + ylab("PCA 3 (2.16%)") +
  xlab("PCA 1 (11.29%)") + geom_point(size=4) + aes(colour = Clusters, shape = Country, stroke = 1) +
  labs(shape = "Country of Origin", colour="Clusters") + scale_color_manual(values = safe_pal_5) +
  scale_shape_manual(values=1:10) + geom_label_repel(aes(label=pca_0.5_merge_K2$'Tag 3'), max.overlaps = Inf)
p
dev.off()

#Shape countries
pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5/PCAs-clusters over 0.7/PCA_0.5_bycountryandcluster_K6_1V3_AM-match-admx-Col.pdf", height=6, width = 6)
p2 <- ggplot(K6_clustering,  aes(x= EV1 , y= EV3)) +
  geom_point() + theme_classic() + ylab("PCA 3 (2.16%)") +
  xlab("PCA 1 (11.29%)") + geom_point(size=4) + aes(colour = Admixed_K6_C_K2, shape = Country, stroke = 1) +
  labs(shape = "Country of Origin", colour="Clusters") + scale_color_manual(values = safe_pal_8) +
  scale_shape_manual(values=1:10) #+ geom_label_repel(aes(label=pca_0.5_merge_K6$'Tag 2'), max.overlaps = Inf)
p2
dev.off()

#Shape photoperiod sensitivity
p2 <- ggplot(K6_clustering,  aes(x= EV1 , y= EV3)) +
  geom_point() + theme_classic() + ylab("PCA 3 (2.16%)") +
  xlab("PCA 1 (11.29%)") + geom_point(size=3) + aes(colour = Admixed_K6_C_K2, shape = `photo sens`, stroke = 1) +
  labs(shape = "Photoperiod sensitivity", colour="Ancestry") + scale_color_manual(values = safe_pal_8)
  #+ geom_label_repel(aes(label=pca_0.5_merge_K6$'Tag 2'), max.overlaps = Inf)
p2

#Shape growth habit
p2 <- ggplot(K6_clustering,  aes(x= EV1 , y= EV3)) +
  geom_point() + theme_classic() + ylab("PCA 3 (2.16%)") +
  xlab("PCA 1 (11.29%)") + geom_point(size=3) + aes(colour = Admixed_K6_C_K2, shape = `Growth habit CIAT then 2023 GH`, stroke = 1) +
  labs(shape = "Growth habit", colour="Ancestry") + scale_color_manual(values = safe_pal_8)
#+ geom_label_repel(aes(label=pca_0.5_merge_K6$'Tag 2'), max.overlaps = Inf)
p2

#K6
p2 <- ggplot(K6_clustering,  aes(x= EV1 , y= EV3)) +
  geom_point() + theme_classic() + ylab("PCA 3 (2.16%)") +
  xlab("PCA 1 (11.29%)") + geom_point(size=3) + aes(colour = Admixed_K6_C_K2) +
  labs(colour="Ancestry") + scale_color_manual(values = safe_pal_8)
#+ geom_label_repel(aes(label=pca_0.5_merge_K6$'Tag 2'), max.overlaps = Inf)
p2

#Merge
pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5/PCAs-clusters over 0.7/PCA_0.5_bycluster_K6_1V2-and-1V3_AM-match-admx-Col.pdf", height=6, width = 10, onefile = F)
p3 <- ggarrange(p1, p2, ncol = 2, nrow = 1, labels = "AUTO", common.legend =  , legend = "right", align = "hv", widths = 4)
p3
dev.off()
?ggarrange

x <- K6_clustering[, c(33,39,43)]

#------------------------------------------------------------------------------------------------------------------
#####Check how PCAs / counrtries look only using CIAT samples - IPK samples may be from markets etcpca_0.5_merge <- merge(pca_0.5, Seed_phenotypes,by = "Sample_name")
pca_0.5_merge_CIAT <- merge(pca_0.5, Seed_phenotypes_CIAT,by = "Sample_name")

pca_0.5_merge_K4_CIAT <- merge(pca_0.5_merge_CIAT, K4_clustering, by = "Sample_name")
pca_0.5_merge_K6_CIAT <- merge(pca_0.5_merge_CIAT, K6_clustering, by = "Sample_name")

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5/PCAs-clusters over 0.7/PCA_0.5_bycountryandcluster_K4_CIAT.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_merge_K4_CIAT,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (2.29%)") +
  xlab("PCA 1 (11.29%)") + geom_point(size=4) + aes(colour = Clusters, shape = Country) +
  labs(colour="Clusters", shape ="Country of Origin") + scale_color_manual(values = safe_pal) + scale_shape_manual(values=1:9)
p
dev.off()

#-------------------------------------------------------------------------------------------------------------------
####### LD 0.5 w.10 - thinned, window size 10  ################
# variance proportion (%)
#pc.percent_0.5_w10 <- pca_0.5_w10$varprop*100
#head(round(pc.percent_0.5_w10, 2))
#21.19  4.15  3.12  2.44  2.03  1.78


pca_0.5_w10 <- 
  read.csv("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Bean project/STATS/22.01 Variant calling - GATK/22.04-filtered_GATK/22.06-Linkage_disequilibrium-in-HPC/PCA_0.5_thin_w10/pca_0.5_w10.txt",
           header=TRUE)

pca_0.5_w10_merge <- merge(pca_0.5_w10, Seed_phenotypes,by = "Sample_name")
pca_0.5_w10_merge_K2 <- merge(pca_0.5_w10_merge, K2_clustering, by = "Sample_name")
pca_0.5_w10_merge_K3 <- merge(pca_0.5_w10_merge, K3_clustering, by = "Sample_name")
pca_0.5_w10_merge_K4 <- merge(pca_0.5_w10_merge, K4_clustering, by = "Sample_name")
pca_0.5_w10_merge_K6 <- merge(pca_0.5_w10_merge, K6_clustering, by = "Sample_name")

#pca_0.5_w10_merge_K6 <- read.csv("../22.11-ADMIXTURE-LD0.5/Q_files/K6_Q_files/K6_clustering_PCA-phenotypes over 0.7.csv")

write.csv(pca_0.5_w10_merge_K6, "pca_0.5_w10_merge_K6")


pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5_thin_w10/PCAs-clusters over 0.7/PCA_0.5_w10_cluster-shape_K6.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_w10_merge_K6,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (4.15%)") +
  xlab("PCA 1 (21.19%)") + geom_point(size=3) + aes(shape = Clusters, stroke = 0.5) +
  labs(shape = "Clusters") +  scale_shape_manual(values=1:7)
p
dev.off()

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5_thin_w10/PCAs-clusters over 0.7/PCA_0.5_w10_cluster-colour_K6.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_w10_merge_K6,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (4.15%)") +
  xlab("PCA 1 (21.19%)") + geom_point(size=3) + aes(colour = Clusters) +
  labs(colour="Clusters") + scale_color_manual(values = safe_pal_7)
p
dev.off()

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5_thin_w10/PCAs-clusters > 0.7/PCA_0.5_w10_bycountryandcluster_K6.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_w10_merge_K6,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (4.15%)") +
  xlab("PCA 1 (21.19%)") + geom_point(size=3) + aes(colour = Country, shape = Clusters) +
  labs(colour="Country of Origin", shape = "Clusters") + scale_color_viridis(discrete=TRUE) +
  scale_shape_manual(values=1:7)
p
dev.off()

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5_thin_w10/PCAs-clusters > 0.7/PCA_0.5_w10_byTypeandcluster_K6.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_w10_merge_K6,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (4.15%)") +
  xlab("PCA 1 (21.19%)") + geom_point(size=3) + aes(colour = Type, shape = Clusters) +
  labs(colour="Type", shape = "Clusters") + scale_color_viridis(discrete=TRUE) +
  scale_shape_manual(values=1:7)
p
dev.off()

#PCA with point labels
pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5_thin_w10/PCAs-clusters > 0.7/PCA_0.5_w10_byCountryandcluster_K6_tags.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_w10_merge_K6,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (4.15%)") +
  xlab("PCA 1 (21.19%)") + geom_point(size=3) + aes(colour = Country, shape = Clusters) +
  labs(colour="Country of Origin", shape = "Clusters") + scale_color_viridis(discrete=TRUE) +
  scale_shape_manual(values=1:7) + geom_label_repel(aes(label=pca_0.5_w10_merge_K6$'Tag 2'), max.overlaps = Inf)
p
dev.off()

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5_thin_w10/PCAs-clusters > 0.7/PCA_0.5_w10_bycountryandcluster_K3_tags.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_w10_merge_K3,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (4.15%)") +
  xlab("PCA 1 (21.19%)") + geom_point(size=3) + aes(colour = Country, shape = Clusters) +
  labs(colour="Country of Origin", shape = "Clusters") + scale_color_viridis(discrete=TRUE) +
  geom_label_repel(aes(label=pca_0.5_w10_merge_K3$'Tag 2'), max.overlaps = Inf)
p
dev.off()

#---------------------------------------------------------
###### Swap shapes and colours ###########

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5_thin_w10/PCAs-clusters over 0.7/PCA_0.5_bycountryandcluster_K2_tags.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_w10_merge_K2,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (4.15%)") +
  xlab("PCA 1 (21.19%)") + geom_point(size=4) + aes(colour = Clusters, shape = Country, stroke = 1) +
  labs(shape = "Country of Origin", colour="Clusters") + scale_color_manual(values = safe_pal_5) +
  scale_shape_manual(values=1:10) + geom_label_repel(aes(label=pca_0.5_w10_merge_K2$'Tag 3'), max.overlaps = Inf)
p
dev.off()

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5_thin_w10/PCAs-clusters over 0.7/PCA_0.5_bycountryandcluster_K6_admxK2.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_w10_merge_K6,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (4.15%)") +
  xlab("PCA 1 (21.19%)") + geom_point(size=4) + aes(colour = Clusters_K2, shape = Country, stroke = 1) +
  labs(shape = "Country of Origin", colour="Clusters") + scale_color_manual(values = safe_pal_8) +
  scale_shape_manual(values=1:10) #+ geom_label_repel(aes(label=pca_0.5_w10_merge_K2$'Tag 2'), max.overlaps = Inf)
p
dev.off()

#---------------------------------------------------------
#Plot 1 v 3 

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5_thin_w10/PCAs-clusters over 0.7/PCA_0.5_bycountryandcluster_K2_1V3.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_w10_merge_K2,  aes(x= EV1 , y= EV3)) +
  geom_point() + theme_classic() + ylab("PCA 3 (3.12%)") +
  xlab("PCA 1 (21.19%)") + geom_point(size=4) + aes(colour = Clusters, shape = Country, stroke = 1) +
  labs(shape = "Country of Origin", colour="Clusters") + scale_color_manual(values = safe_pal_5) +
  scale_shape_manual(values=1:10) #+ geom_label_repel(aes(label=pca_0.5_w10_merge_K4$'Tag 2'), max.overlaps = Inf)
p
dev.off()

pdf("22.06-Linkage_disequilibrium-in-HPC/PCA_0.5_thin_w10/PCAs-clusters over 0.7/PCA_0.5_bycountryandcluster_K6_1V3.pdf", height=6, width = 6)
p <- ggplot(pca_0.5_w10_merge_K6,  aes(x= EV1 , y= EV3)) +
  geom_point() + theme_classic() + ylab("PCA 3 (3.12%)") +
  xlab("PCA 1 (21.19%)") + geom_point(size=4) + aes(colour = Clusters, shape = Country, stroke = 1) +
  labs(shape = "Country of Origin", colour="Clusters") + scale_color_manual(values = safe_pal_7) +
  scale_shape_manual(values=1:10) #+ geom_label_repel(aes(label=pca_0.5_w10_merge_K6$'Tag 2'), max.overlaps = Inf)
p
dev.off()

p1 <- 

#------------------------------------------------------------------------------------------------------------------
####### Clusters for K2 #################
#Select Individuals or sample names from K2 Clusters 
#I will include Cluster 1 and admixed then Cluster 2 and admixed 
#The idea is to create 2 files for each cluster. I will use this for PCAs and to run admixture again on each cluster. 
#Cluster_join <- join(samples_uniq_all, K2_clustering)
library(readr)
samples_uniq_all <- read_csv("samples_uniq_all.txt", col_names = FALSE)
samples_uniq_all <- as.data.frame(samples_uniq_all)

K2_clustering$sample_uniq_all <- samples_uniq_all

Cluster_1 = dplyr::filter(K2_clustering, Clusters %in% c('Cluster 1', 'Admixed'))
Cluster_2 = dplyr::filter(K2_clustering, Clusters %in% c('Cluster 2', 'Admixed'))

Cluster_1_names = as_data_frame(Cluster_1$sample_uniq_all)
Cluster_2_names = as_data_frame(Cluster_2$sample_uniq_all)

write.csv(
  file = "Bean project/STATS/22.01 Variant calling - GATK/22.11-ADMIXTURE-LD0.5/Q_files/K2_Q_files/K2_Cluster_2_sample_uniq_names.csv", 
  Cluster_2_names
)


###PCAs-------------------------------------------------------------------
# variance proportion (%)
#pc.percent <- pca__ADMX_K2_C1$varprop*100
#7.57 5.36 3.12 3.05 2.98 2.73

# variance proportion (%)
#pc.percent <- pca__ADMX_K2_C2$varprop*100
#10.09  6.99  4.85  3.94  3.09  2.11

K2_C1_clustering <- 
  read.csv(
    "Bean project/STATS/22.01 Variant calling - GATK/22.04-filtered_GATK/22.12-ADMX-K=2-vcfs/Q_files/C1/K2_Q_files/K2_Clusters over 0.7.csv"
  )

K5_C2_clustering <- 
  read.csv(
    "/Bean project/STATS/22.01 Variant calling - GATK/22.04-filtered_GATK/22.12-ADMX-K=2-vcfs/Q_files/C2/K5_Q_files/K5_Clusters over 0.7.csv"
  )

K7_C2_clustering <- 
  read.csv(
    "Bean project/STATS/22.01 Variant calling - GATK/22.04-filtered_GATK/22.12-ADMX-K=2-vcfs/Q_files/C2/K7_Q_files/K7_Clusters over 0.7.csv"
  )

names(K7_C2_clustering)[names(K7_C2_clustering) == 'Sample.name'] <- 'Sample_name'

pca_K2_C1 <- 
  read.csv(
    "Bean project/STATS/22.01 Variant calling - GATK/22.04-filtered_GATK/22.12-ADMX-K=2-vcfs/pca_K2_C1.txt",
           header=TRUE)

pca_K2_C2 <- 
  read.csv(
    "~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Bean project/STATS/22.01 Variant calling - GATK/22.04-filtered_GATK/22.12-ADMX-K=2-vcfs/pca_K2_C2.txt",
    header=TRUE)

pca_K2_C2_merge <- merge(pca_K2_C2, Seed_phenotypes, by = "Sample_name")
pca_K2_C2_merge_K7 <- merge(pca_K2_C2_merge, K7_C2_clustering, by = "Sample_name")

write.csv(
  file = "Bean project/STATS/22.01 Variant calling - GATK/22.04-filtered_GATK/22.12-ADMX-K=2-vcfs/Q_files/C2/K5_Q_files/K2_C2_clustering_K7_PCA-phenotypes > 0.7.csv", 
  pca_K2_C2_merge_K7 
)

###Plots 
pdf("Bean project/STATS/22.01 Variant calling - GATK/22.04-filtered_GATK/22.12-ADMX-K=2-vcfs/Q_files/C1/PCAs/PCA_K2_C1_bycountryandcluster_K2.pdf", height=6, width = 6)
p <- ggplot(pca_K2_C1_merge_K2,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (5.36%)") +
  xlab("PCA 1 (7.57%)") + geom_point(size=4) + aes(colour = Clusters, shape = Country, stroke = 1) +
  labs(shape = "Country of Origin", colour="Clusters") + scale_color_manual(values = safe_pal_5) +
  scale_shape_manual(values=1:10) #+ geom_label_repel(aes(label=pca_K2_C1_merge_K2$'Tag 3'), max.overlaps = Inf)
p
dev.off()

pdf("Bean project/STATS/22.01 Variant calling - GATK/22.04-filtered_GATK/22.12-ADMX-K=2-vcfs/Q_files/C1/PCAs/PCA_K2_C1_bycountryandcluster_K2_1V3.pdf", height=6, width = 6)
p <- ggplot(pca_K2_C1_merge_K2,  aes(x= EV1 , y= EV3)) +
  geom_point() + theme_classic() + ylab("PCA 3 (3.12%)") +
  xlab("PCA 1 (7.57%)") + geom_point(size=4) + aes(colour = Clusters, shape = Country, stroke = 1) +
  labs(shape = "Country of Origin", colour="Clusters") + scale_color_manual(values = safe_pal_5) +
  scale_shape_manual(values=1:10) #+ geom_label_repel(aes(label=pca_K2_C1_merge_K2$'Tag 3'), max.overlaps = Inf)
p
dev.off()

##C2
pdf("Bean project/STATS/22.01 Variant calling - GATK/22.04-filtered_GATK/22.12-ADMX-K=2-vcfs/Q_files/C2/PCAs/PCA_K2_C2_bycountryandcluster_K7.pdf", height=6, width = 6)
p <- ggplot(pca_K2_C2_merge_K7,  aes(x= EV1 , y= EV2)) +
  geom_point() + theme_classic() + ylab("PCA 2 (6.99%)") +
  xlab("PCA 1 (10.09%)") + geom_point(size=4) + aes(colour = Clusters, shape = Country, stroke = 1) +
  labs(shape = "Country of Origin", colour="Clusters") + scale_color_manual(values = safe_pal_8) +
  scale_shape_manual(values=1:10) #+ geom_label_repel(aes(label=pca_K2_C2_merge_K7$'Tag 3'), max.overlaps = Inf)
p
dev.off()

pdf("Bean project/STATS/22.01 Variant calling - GATK/22.04-filtered_GATK/22.12-ADMX-K=2-vcfs/Q_files/C2/PCAs/PCA_K2_C1_bycountryandcluster_K2_1V3.pdf", height=6, width = 6)
p <- ggplot(pca_K5_C2_merge_K5,  aes(x= EV1 , y= EV3)) +
  geom_point() + theme_classic() + ylab("PCA 3 (4.85%)") +
  xlab("PCA 1 (10.09%)") + geom_point(size=4) + aes(colour = Clusters, shape = Country, stroke = 1) +
  labs(shape = "Country of Origin", colour="Clusters") + scale_color_manual(values = safe_pal_7) +
  scale_shape_manual(values=1:10) #+ geom_label_repel(aes(label=pca_K5_C2_merge_K5$'Tag 3'), max.overlaps = Inf)
p
dev.off()

------------------------------------------------------------------------------------------------------------------
###### Check whether samples in Mesoamerican clusters match 
#This is K=4 Cluster 4 and K=6 Cluster 3 and5 

C3_K4 <- dplyr::filter(K4_clustering, Clusters %in% c('Cluster 4'))
C35_K6 <- dplyr::filter(K6_clustering, Clusters %in% c('Cluster 3', 'Cluster 5'))

K4 <- as_data_frame(C3_K4$Sample_name)
K6 <- as_data_frame(C35_K6$Sample_name)

K4_K6 <- merge(C35_K6,C3_K4, by ='Sample_name')

install.packages('sqldf')
library('sqldf')
require(sqldf)
res <- sqldf('SELECT * FROM K4 EXCEPT SELECT * FROM K6')
