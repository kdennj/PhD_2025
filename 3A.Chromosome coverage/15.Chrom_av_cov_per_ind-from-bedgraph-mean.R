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

rm(list = ls())

setwd("STATS/22.03 Chromosome coverage/bowtie/23.02-bedgraphs-mean/")

data <- read.csv("S0001_G50516H_chr_cov.csv")
str(data)
dim(data)

data <- data[1:22, ]
chroms <- str_pad(1:11, pad = 0,width = 2 , "left")

##Add column for chromosome numbers
data$Chromosome <- chroms

#Add columns for reference genome
Andean <- data[data$chrom %like% "ChrA",]
Andeans <- as.list(Andean$X)
data$genome <- "M"
data$genome[data$X %in% Andeans] = "A"


mycolors_change <- c("dodgerblue","firebrick3")

###Plot for one individual - testing - files from bedgraph mean 
p <- ggplot(data, aes(Chromosome, cov))
p1 <- p + geom_line(col = "grey60")
p2 <- p1 + geom_point(aes(colour=genome))
p3 <- p2 + theme(axis.text.x=element_text(angle = 90, hjust = 1, size=6)) + xlab("Chromosome") + ylab("mean percentage coverage mapped to each chromosome")
p4 <- p3 +  scale_color_manual(values=mycolors_change)
p4

pdf(file = "plots_ind_av_chrom/S0001_G50516H_chr_cov_av.pdf",width=6,height=6)
p4
dev.off()

###For loop
pop <- read.table("../bowtie_practise.txt") #list of sample names
slist <- pop$V1
mycolors_change <- c("dodgerblue","firebrick3")

for (sam in slist){
  name <- paste(sam,"_chr_cov.csv", sep = "")
  print (name) 
  
  data <- read.csv(name, header = TRUE)
  data <- data[1:22, ]
  
  ##Add column for chromosome numbers
  chroms <- str_pad(1:11, pad = 0,width = 2 , "left")
  data$Chromosome <- chroms
  
  #Add columns for reference genome
  Andean <- data[data$chrom %like% "ChrA",]
  Andeans <- as.list(Andean$X)
  data$genome <- "M"
  data$genome[data$X %in% Andeans] = "A"
  
  ###Plot for one individual - testing - files from bedgraph mean 
  aname <- paste(sam,"_chr_cov_av.pdf",sep = "")
  print(aname)
  
  p <- ggplot(data, aes(Chromosome, cov))
  p1 <- p + geom_line(col = "grey60")
  p2 <- p1 + geom_point(aes(colour=genome))
  p3 <- p2 + theme(axis.text.x=element_text(angle = 90, hjust = 1, size=6)) + xlab("Chromosome") + ylab("mean percentage coverage mapped to each chromosome")
  p4 <- p3 +  scale_color_manual(values=mycolors_change)
  p4
  
  pdf(file = paste("plots_ind_av_chrom/", aname, sep = ""), width =7, height = 7)
  plot(p4)
  dev.off()
  
}

#####Plot based on subgroups 
###Need 1 file with same information but for all samples

##Add column for sample names
data$samples <- "S0001_G50516H"
data[c('Number', 'Sample_name')] <- str_split_fixed(data$samples, "_", 2)

data_all <- data.frame()

###For loop
pop <- read.table("../bowtie_practise.txt") #list of sample names
slist <- pop$V1
mycolors_change <- c("dodgerblue","firebrick3")


for (sam in slist){
  name <- paste(sam,"_chr_cov.csv", sep = "")
  print (name) 
  
  data <- read.csv(name, header = TRUE)
  data <- data[1:22, ]
  
  #Add columns with the names
  data$samples <- sam
  data[c('Number', 'Sample_name')] <- str_split_fixed(data$samples, "_", 2)
  
  ##Add column for chromosome numbers
  chroms <- str_pad(1:11, pad = 0,width = 2 , "left")
  data$Chromosome <- chroms
  
  #Add columns for reference genome
  Andean <- data[data$chrom %like% "ChrA",]
  Andeans <- as.list(Andean$X)
  data$genome <- "M"
  data$genome[data$X %in% Andeans] = "A"
  
  
  if (sam == "S0001_G50516H"){
  data_all <- data
  }
  
  else {
   data_all <- rbind(data, data_all)
  }
}

###Plot for many individuals - testing - files from bedgraph mean 
p <- ggplot(data_all, aes(Sample_name, cov))
p1 <- p + geom_line(col = "grey60")
p2 <- p1 + geom_point(aes(colour=genome))
p3 <- p2 + theme(axis.text.x=element_text(angle = 90, hjust = 1, size=6)) + xlab("Chromosome") + ylab("mean percentage coverage mapped to each chromosome")
p4 <- p3 +  scale_color_manual(values=mycolors_change)
p5 <- p4 + facet_wrap(~Chromosome)
p5

pdf(file = "plots_av_chrom_grouped/practise_3_chr_cov_av.pdf",width=8,height=7)
p5
dev.off()


#############Repeat but for bwamem files------------------------------------------------
setwd("~/OneDrive - Norwich BioScience Institutes/Bean project/STATS/22.03 Chromosome coverage/bwa/23.01-mean_cov")

data <- read.csv("S0001_G50516H_chr_cov.csv")
str(data)
dim(data)

data <- data[1:22, ]
chroms <- str_pad(1:11, pad = 0,width = 2 , "left")

##Add column for chromosome numbers
data$Chromosome <- chroms

#Add columns for reference genome
Andean <- data[data$chrom %like% "ChrA",]
Andeans <- as.list(Andean$X)
data$genome <- "M"
data$genome[data$X %in% Andeans] = "A"


mycolors_change <- c("dodgerblue","firebrick3")

###Plot for one individual - testing - files from bedgraph mean 
p <- ggplot(data, aes(Chromosome, cov))
p1 <- p + geom_line(col = "grey60")
p2 <- p1 + geom_point(aes(colour=genome))
p3 <- p2 + theme(axis.text.x=element_text(angle = 90, hjust = 1, size=6)) + xlab("Chromosome") + ylab("mean percentage coverage mapped to each chromosome")
p4 <- p3 +  scale_color_manual(values=mycolors_change)
p4

pdf(file = "plots_ind_av_chrom/S0001_G50516H_chr_cov_av.pdf",width=6,height=6)
p4
dev.off()

###For loop

pop <- read.table("../../samples_uniq_all.txt") #list of sample names
slist <- pop$V1
mycolors_change <- c("dodgerblue","firebrick3")

for (sam in slist){
  name <- paste(sam,"_chr_cov.csv", sep = "")
  print (name) 
  
  data <- read.csv(name, header = TRUE)
  data <- data[1:22, ]
  
  ##Add column for chromosome numbers
  chroms <- str_pad(1:11, pad = 0,width = 2 , "left")
  data$Chromosome <- chroms
  
  #Add columns for reference genome
  Andean <- data[data$chrom %like% "ChrA",]
  Andeans <- as.list(Andean$X)
  data$genome <- "M"
  data$genome[data$X %in% Andeans] = "A"
  
  ###Plot for one individual - testing - files from bedgraph mean 
  aname <- paste(sam,"_chr_cov_av.pdf",sep = "")
  print(aname)
  
  p <- ggplot(data, aes(Chromosome, cov))
  p1 <- p + geom_line(col = "grey60")
  p2 <- p1 + geom_point(aes(colour=genome))
  p3 <- p2 + theme(axis.text.x=element_text(angle = 90, hjust = 1, size=6)) + xlab("Chromosome") + ylab("mean percentage coverage mapped to each chromosome")
  p4 <- p3 +  scale_color_manual(values=mycolors_change)
  p4
  
  pdf(file = paste("plots_ind_av_chrom/", aname, sep = ""), width =7, height = 7)
  plot(p4)
  dev.off()
  
}






