#Install packages
#install.packages("ggplot2")
library(ggplot2)
library(dplyr)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("viridis")
library(viridis)
#install.packages("hrbrthemes")
#library(hrbrthemes)

library("gridExtra")
library("cowplot")
library(patchwork)
#install.packages("ggpubr")
library("ggpubr")
library("data.table")

#Import data
rm(list = ls())

setwd("/STATS/22.01 Variant calling - GATK/22.04-filtered_GATK")

library(readxl)
Summary <- read_excel("/STATS/22.01 Variant calling - GATK/22.04-filtered_GATK/Taxa Summary of All_Andean_Haplotypecaller_filtered.xlsx")
#View(Summary)

names(Summary)[names(Summary) == 'Taxa Name'] <- 'TAXA'
Summary_new <- separate(Summary, TAXA, sep = "_", into = c("No", "Sample_name"))
Summary_new$Sample_name[45] <- 'G50516I' 

Seed_phenotypes<- read_excel("Common bean samples/22.02 Seed phenotypes/22.02-updated-phenotypes-R.xlsx")
#View(Seed_phenotypes)

JOIN <- merge(x=Summary_new, y=Seed_phenotypes, by="Sample_name")

safe_pal_4 <- c("#1b9e77", "#d95f02", "#7570b3", "#000000")
safe_pal_8 <- c("#000000", "#888888", "#5EC962", "#117733", "#332288", "#6699CC", "#CC5500","#AA4499")
safe_pal_8 <- c("#117733", "#000000", "#888888", "#332288", "#5EC962", "#AA4499", "#6699CC", "#CC5500")

#Mesoamerican
#Colour
M <- ggplot(JOIN, aes(x= `%M`, y= `Proportion Heterozygous`)) +
  geom_point() + theme_classic() + xlab("% Mesoamerican alignment ") +
  ylab("Proportion of heterozygous sites") + geom_point(size=2.7) + aes(colour = `Major and second`, shape = Type) +
  labs(colour="Major Seed Colour", shape="Type") #+ scale_color_manual(values = safe_pal_8)  + 
  theme(legend.justification = c(1, 0),
        legend.position = c(1,0),) 
M

pdf(file = "Meso_hetero_accession_major_seed_colour.pdf",width = 7, height = 6)
plot(M)
dev.off()


#Andean
A <- ggplot(JOIN, aes(x= `%A`, y= `Proportion Heterozygous`, colour = Admixed_K6_C_K2)) +
  geom_point() + theme_classic() + xlab("% Andean alignment ") +
  ylab("Proportion of heterozygous sites") + geom_point(size=2.5) +
  labs(colour="Ancestry") + 
  scale_color_manual(values = safe_pal_8)
  #scale_color_brewer(palette = "Dark2")
  

#theme(legend.justification = c(1, 0),
        #legend.position = c(1,0),) 
A

pdf(file = "Andean_hetero_accession_K6_ADMX.pdf",width = 7, height = 6)
plot(A)
dev.off()



##Arrange
p3 = grid.arrange(A, M + theme(legend.position="none"))

p5 = ggarrange(A, M, ncol = 2, nrow = 1, legend = "bottom", common.legend = T)
p5

?ggarrange

pdf(file = "A_M_hetero_accession_Liborino.pdf",width = 10, height = 5)
plot(p5)
dev.off()
