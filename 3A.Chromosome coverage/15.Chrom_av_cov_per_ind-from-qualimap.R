library(reshape)
library(ggplot2)
library(stringr)
library("plyr")
library(dplyr)
library(tidyr)
library(tidyverse)
library(remotes)
library("devtools")
library("tidyselect")
library("data.table")
library(readxl)

rm(list = ls())

setwd("STATS/22.03 Chromosome coverage/bowtie2/23.02-qualimap-res/")


data <- read.delim("S0001_G50516H_genome_results_tail.txt", header=FALSE)
str(data)
dim(data)

data <- data[,2:6]

##Change column names
colnames(data) <- c("Chrom", "Length", "Mapped_bases", "Mean_cov", "SD")

##Remove scaffold columns
data <- data[data$Chrom %like% "Chr",]

##Add column for chromosome numbers
chroms <- str_pad(1:11, pad = 0,width = 2 , "left")
data$Chromosome <- chroms


#Add columns for reference genome
X <- c(1:22)
data$X <- X
Andean <- data[data$Chrom %like% "ChrA",]
Andeans <- as.list(Andean$X)
data$genome <- "M"
data$genome[data$X %in% Andeans] = "A"

mycolors_change <- c("dodgerblue","firebrick3")

###Plot for one individual - testing - files from bedgraph mean 
p <- ggplot(data, aes(Chromosome, Mean_cov))
p1 <- p + geom_line(col = "grey60")
p2 <- p1 + geom_point(aes(colour=genome))
p3 <- p2 + theme(axis.text.x=element_text(angle = 90, hjust = 1, size=6)) + xlab("Chromosome") + ylab("Mean coverage mapped to each chromosome from qualimap")
p4 <- p3 +  scale_color_manual(values=mycolors_change)
p4

pdf(file = "plots_ind_av_mean/S0001_G50516H_chr_cov_av.pdf",width=6,height=6)
p4
dev.off()

###For loop --------------- --------------- --------------- ---------------
pop <- read.table("../../samples_uniq_all.txt") #list of sample names
slist <- pop$V1
mycolors_change <- c("dodgerblue","firebrick3")

for (sam in slist){
  name <- paste(sam,"_genome_results_tail.txt", sep = "")
  print (name) 
  
  data <- read.delim(name, header = FALSE)
  
  data <- data[,2:6]
  
  ##Change column names
  colnames(data) <- c("Chrom", "Length", "Mapped_bases", "Mean_cov", "SD")
  
  ##Remove scaffold columns
  data <- data[data$Chrom %like% "Chr",]
  
  ##Add column for chromosome numbers
  chroms <- str_pad(1:11, pad = 0,width = 2 , "left")
  data$Chromosome <- chroms
  
  #Add columns for reference genome
  X <- c(1:22)
  data$X <- X
  Andean <- data[data$Chrom %like% "ChrA",]
  Andeans <- as.list(Andean$X)
  data$genome <- "M"
  data$genome[data$X %in% Andeans] = "A"

  ###Plot for one individual - testing 
  aname <- paste(sam,"_chr_cov_av.pdf",sep = "")
  print(aname)
  
  p <- ggplot(data, aes(Chromosome, Mean_cov))
  p1 <- p + geom_line(col = "grey60")
  p2 <- p1 + geom_point(aes(colour=genome))
  p3 <- p2 + theme(axis.text.x=element_text(angle = 90, hjust = 1, size=6)) + xlab("Chromosome") + ylab("Mean coverage mapped to each chromosome from qualimap")
  p4 <- p3 +  scale_color_manual(values=mycolors_change)
  p4
  
  pdf(file = paste("plots_ind_av_mean/", aname, sep = ""), width =7, height = 7)
  plot(p4)
  dev.off()

}

#####Plot based on subgroups 
###Need 1 file with same information but for all samples

###For loop
pop <- read.table("../../bowtie_practise.txt") #list of sample names
slist <- pop$V1
mycolors_change <- c("dodgerblue","firebrick3")


for (sam in slist){
  name <- paste(sam,"_genome_results_tail.txt", sep = "")
  print (name) 
  
  data <- read.delim(name, header = FALSE)
  
  data <- data[,2:6]
  
  ##Change column names
  colnames(data) <- c("Chrom", "Length", "Mapped_bases", "Mean_cov", "SD")
  
  ##Remove scaffold columns
  data <- data[data$Chrom %like% "Chr",]
  
  ##Add column for chromosome numbers
  chroms <- str_pad(1:11, pad = 0,width = 2 , "left")
  data$Chromosome <- chroms
  
  #Add columns for reference genome
  X <- c(1:22)
  data$X <- X
  Andean <- data[data$Chrom %like% "ChrA",]
  Andeans <- as.list(Andean$X)
  data$genome <- "M"
  data$genome[data$X %in% Andeans] = "A"
  
  #Add columns with the names
  data$samples <- sam
  data[c('Number', 'Sample_name')] <- str_split_fixed(data$samples, "_", 2)
  
  if (sam == "S0001_G50516H"){
    data_all <- data
  }
  
  else {
    data_all <- rbind(data, data_all)
  }
}

data_all <- data_all %>% arrange(samples)

write.csv(data_all, file = "bowtie2_all.csv")

data_all <- read.csv("bowtie2_all.csv")

###Plot for many individuals - testing - files from bedgraph mean 
p <- ggplot(data_all, aes(Sample_name, Mean_cov))
p1 <- p + geom_line(col = "grey60")
p2 <- p1 + geom_point(aes(colour=genome))
p3 <- p2 + theme(axis.text.x=element_text(angle = 90, hjust = 1, size=6)) + xlab("Chromosome") + ylab("mean percentage coverage mapped to each chromosome")
p4 <- p3 +  scale_color_manual(values=mycolors_change)
p5 <- p4 + facet_wrap(~Chromosome)
p5

pdf(file = "plots_av_chrom_grouped/practise_3_chr_cov_av.pdf",width=9,height=9)
plot(p5)
dev.off()

###Grouping y ADMX K=6 -----------------------------
#####Rename G505161 TO G50516I
data_all$samples[969:990] <- 'S0046_G50516I'
data_all$Sample_name[969:990] <- 'G50516I'

##### Based on K6 Clustering_admixture for bwa

K6_clustering <- 
  read_excel("STATS/22.01 Variant calling - GATK/22.11-ADMIXTURE-LD0.5/Q_files/K6_Q_files/K6_clustering_PCA-phenotypes over 0.7 ADMX.xlsx")

K6_clustering$ClustersADMX <- str_replace_all(K6_clustering$ClustersADMX, "Cluster", "Ancestry")

K6_clustering <- K6_clustering[,c(2, 33, 10)]

K6_data_all <- merge(K6_clustering, data_all, by ="Sample_name")

##Split into clusters
Cluster1 <- K6_data_all %>% filter(grepl('Cluster 1', ClustersADMX))
Cluster2 <- K6_data_all %>% filter(grepl('Cluster 2', ClustersADMX))
Cluster3 <- K6_data_all %>% filter(grepl('Cluster 3', ClustersADMX))
Cluster4 <- K6_data_all %>% filter(grepl('Cluster 4', ClustersADMX))
Cluster5 <- K6_data_all %>% filter(grepl('Cluster 5', ClustersADMX))
Cluster6 <- K6_data_all %>% filter(grepl('Cluster 6', ClustersADMX))
ADMIXED_A <- K6_data_all %>% filter(grepl('Admixed_A', ClustersADMX))
ADMIXED_M <- K6_data_all %>% filter(grepl('Admixed_M', ClustersADMX))
ADMIXED_AM <- K6_data_all %>% filter(grepl('Admixed-AM', ClustersADMX))

Colombia <- K6_data_all %>% filter(grepl('Colombia', Country))
Colombia <- unique(Colombia$Sample_name)

###Plot clusters
Cluster <- Cluster6
Cluster <- ADMIXED_AM

p <- ggplot(Cluster, aes(Sample_name, Mean_cov))
p1 <- p + geom_line(col = "grey60")
p2 <- p1 + geom_point(aes(colour=genome))
p3 <- p2 + theme(axis.text.x=element_text(angle = 90, hjust = 1, size=6)) + xlab("Chromosome") + ylab("mean percentage coverage mapped to each chromosome")
p4 <- p3 +  scale_color_manual(values=mycolors_change)
p5 <- p4 + facet_wrap(~Chromosome) 
p5

pdf(file = "plots_av_chrom_grouped/ADMIXED_AM.pdf",width=9,height=8)
pdf(file = "plots_av_chrom_grouped/Cluster6.pdf",width=9,height=8)
p5
dev.off()

####Plot all samples and wrap by chromosome ############
K6_data_all <- K6_data_all %>% arrange(Admixed_K6_C)

Chrom1 <- K6_data_all[K6_data_all$Chromosome == "1",]
Chrom2 <- K6_data_all[K6_data_all$Chromosome == "2",]
Chrom3 <- K6_data_all[K6_data_all$Chromosome == "3",]
Chrom4 <- K6_data_all[K6_data_all$Chromosome == "4",]
Chrom5 <- K6_data_all[K6_data_all$Chromosome == "5",]
Chrom6 <- K6_data_all[K6_data_all$Chromosome == "6",]
Chrom7 <- K6_data_all[K6_data_all$Chromosome == "7",]
Chrom8 <- K6_data_all[K6_data_all$Chromosome == "8",]
Chrom9 <- K6_data_all[K6_data_all$Chromosome == "9",]
Chrom10 <- K6_data_all[K6_data_all$Chromosome == "10",]
Chrom11 <- K6_data_all[K6_data_all$Chromosome == "11",]

Chrom <- Chrom %>% arrange(Admixed_K6_C)
ADMX <- unique(K6_data_all$Sample_name)

Chrom <- Chrom11

print(Colombia)

p <- ggplot(Chrom, aes(x=factor(Sample_name, level = ADMX) , Mean_cov)) 
p1 <- p + geom_line(col = "grey60")
p2 <- p1 + geom_point(aes(colour=genome)) 
p3 <- p2 + theme(axis.text.x=element_text(angle = 90, hjust = 1, size=6)) + xlab("Sample name") + ylab("Mean percentage coverage mapped to each chromosome")
p4 <- p3 +  scale_color_manual(values=mycolors_change)
p5 <- p4 + facet_grid(~Admixed_K6_C, scales = 'free') + ggtitle("Chromosome 11") #+ scale_x_discrete(labels = ~ ifelse(.x == "G11819", paste0(.x, "*"), .x))
p5

pdf(file = "plots_all/New-clusters/Chromosome11_ADMX-Col.pdf",width=25,height=8)
p5
dev.off()


#############Repeat but for bwamem files------------------------------------------------
setwd("/STATS/22.03 Chromosome coverage/bwa/23.02-qualimap-res")
rm(list = ls())

###For loop --------------- --------------- --------------- ---------------
pop <- read.table("../../samples_uniq_all.txt") #list of sample names
slist <- pop$V1
mycolors_change <- c("dodgerblue","firebrick3")

for (sam in slist){
  name <- paste(sam,"_genome_results_tail.txt", sep = "")
  print (name) 
  
  data <- read.delim(name, header = FALSE)
  
  data <- data[,2:6]
  
  ##Change column names
  colnames(data) <- c("Chrom", "Length", "Mapped_bases", "Mean_cov", "SD")
  
  ##Remove scaffold columns
  data <- data[data$Chrom %like% "Chr",]
  
  ##Add column for chromosome numbers
  chroms <- str_pad(1:11, pad = 0,width = 2 , "left")
  data$Chromosome <- chroms
  
  #Add columns for reference genome
  X <- c(1:22)
  data$X <- X
  Andean <- data[data$Chrom %like% "ChrA",]
  Andeans <- as.list(Andean$X)
  data$genome <- "M"
  data$genome[data$X %in% Andeans] = "A"
  
  ###Plot for one individual - testing - files from bedgraph mean 
  aname <- paste(sam,"_chr_cov_av.pdf",sep = "")
  print(aname)
  
  p <- ggplot(data, aes(Chromosome, Mean_cov))
  p1 <- p + geom_line(col = "grey60")
  p2 <- p1 + geom_point(aes(colour=genome))
  p3 <- p2 + theme(axis.text.x=element_text(angle = 90, hjust = 1, size=6)) + xlab("Chromosome") + ylab("Mean coverage mapped to each chromosome from qualimap")
  p4 <- p3 +  scale_color_manual(values=mycolors_change)
  p4
  
  pdf(file = paste("plots_ind_av_mean/", aname, sep = ""), width =7, height = 7)
  plot(p4)
  dev.off()
  
}

##### For loop to create file with sample names -------------------------------
pop <- read.table("../../samples_uniq_all.txt") #list of sample names
slist <- pop$V1
mycolors_change <- c("dodgerblue","firebrick3")

for (sam in slist){
  name <- paste(sam,"_genome_results_tail.txt", sep = "")
  print (name) 
  
  data <- read.delim(name, header = FALSE)
  
  data <- data[,2:6]
  
  ##Change column names
  colnames(data) <- c("Chrom", "Length", "Mapped_bases", "Mean_cov", "SD")
  
  ##Remove scaffold columns
  data <- data[data$Chrom %like% "Chr",]
  
  ##Add column for chromosome numbers
  chroms <- str_pad(1:11, pad = 0,width = 2 , "left")
  data$Chromosome <- chroms
  
  #Add columns for reference genome
  X <- c(1:22)
  data$X <- X
  Andean <- data[data$Chrom %like% "ChrA",]
  Andeans <- as.list(Andean$X)
  data$genome <- "M"
  data$genome[data$X %in% Andeans] = "A"
  
  #Add columns with the names
  data$samples <- sam
  data[c('Number', 'Sample_name')] <- str_split_fixed(data$samples, "_", 2)
  
  if (sam == "S0001_G50516H"){
    data_all <- data
  }
  
  else {
    data_all <- rbind(data, data_all)
  }
}


write.csv(data_all, file = "All_data_genome_files_tail.csv")

data_all <- read.csv("All_data_genome_files_tail.csv", header = TRUE)

#####Rename G505161 TO G50516I
data_all$samples[2179:2200] <- 'S0046_G50516I'
data_all$Sample_name[2179:2200] <- 'G50516I'

##### Based on K6 Clustering_admixture for bwa

K6_clustering <- 
  read.csv("/STATS/22.03 Chromosome coverage/K6_clustering_PCA-phenotypes over 0.7 ADMX.csv",
           header=TRUE, )

K6_clustering$ClustersADMX <- str_replace_all(K6_clustering$ClustersADMX, "Cluster", "Ancestry")

K6_clustering <- K6_clustering[,c(3, 33)]

K6_data_all <- merge(K6_clustering, data_all, by ="Sample_name")

##Split into clusters
Cluster1 <- K6_data_all %>% filter(grepl('Cluster 1', ClustersADMX))
Cluster2 <- K6_data_all %>% filter(grepl('Cluster 2', ClustersADMX))
Cluster3 <- K6_data_all %>% filter(grepl('Cluster 3', ClustersADMX))
Cluster4 <- K6_data_all %>% filter(grepl('Cluster 4', ClustersADMX))
Cluster5 <- K6_data_all %>% filter(grepl('Cluster 5', ClustersADMX))
Cluster6 <- K6_data_all %>% filter(grepl('Cluster 6', ClustersADMX))
ADMIXED_A <- K6_data_all %>% filter(grepl('Admixed_A', ClustersADMX))
ADMIXED_M <- K6_data_all %>% filter(grepl('Admixed_M', ClustersADMX))
ADMIXED_AM <- K6_data_all %>% filter(grepl('Admixed-AM', ClustersADMX))

###Plot clusters
Cluster <- Cluster6
Cluster <- ADMIXED_M

p <- ggplot(Cluster, aes(Sample_name, Mean_cov))
p1 <- p + geom_line(col = "grey60")
p2 <- p1 + geom_point(aes(colour=genome))
p3 <- p2 + theme(axis.text.x=element_text(angle = 90, hjust = 1, size=6)) + xlab("Chromosome") + ylab("mean percentage coverage mapped to each chromosome")
p4 <- p3 +  scale_color_manual(values=mycolors_change)
p5 <- p4 + facet_wrap(~Chromosome)
p5

pdf(file = "plot_clusters_av_mean/Cluster6.pdf",width=9,height=8)
pdf(file = "plot_clusters_av_mean/ADMIXED_M.pdf",width=9,height=8)
p5
dev.off()


