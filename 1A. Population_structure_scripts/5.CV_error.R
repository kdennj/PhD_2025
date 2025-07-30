# install.packages('devtools',dependencies=T)
library(devtools)
# install_github('royfrancis/pophelper')
library(pophelper)

# packageDescription("pophelper", fields="Version")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(tidyr)
library("RColorBrewer")
library(remotes)
# install.packages("reshape2")
library(reshape)
library(reshape2)
library(plyr)
library(dplyr)


## CV error for ADMIXTURE ###############
all_Ks_uniq_log <- read.table("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Bean project/STATS/22.01 Variant calling - GATK/22.11-ADMIXTURE-LD0.5/Log_files/all_Ks_uniq_log.txt", quote = "\"", comment.char = "")

names(all_Ks_uniq_log)[1] <- paste("K")
names(all_Ks_uniq_log)[2] <- paste("Cross-validation error")

plot(all_Ks_uniq_log)

pdf("../Log_files/CV-error-plot.pdf", height = 5, width = 5)
p <- ggplot(all_Ks_uniq_log, aes(x = K, y = `Cross-validation error`)) +
    geom_point() +
    theme_classic() +
    geom_point(size = 2) +
    geom_line()
p
dev.off()

## Extract clusters for PCA/ tree colouring

## CV error on ADMIXTURE on cluster 1 and 2
setwd("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Bean project/STATS/22.01 Variant calling - GATK/22.04-filtered_GATK/22.12-ADMX-K=2-vcfs/")

all_Ks_C1_uniq <- read.table("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Bean project/STATS/22.01 Variant calling - GATK/22.04-filtered_GATK/22.12-ADMX-K=2-vcfs/Log-files/all_Ks_C1_uniq.txt",
    quote = "\"", comment.char = ""
)
all_Ks_C2_uniq <- read.table("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Bean project/STATS/22.01 Variant calling - GATK/22.04-filtered_GATK/22.12-ADMX-K=2-vcfs/Log-files/all_Ks_C2_uniq.txt",
    quote = "\"", comment.char = ""
)

names(all_Ks_C1_uniq)[1] <- paste("K")
names(all_Ks_C1_uniq)[2] <- paste("Cross-validation error")

names(all_Ks_C2_uniq)[1] <- paste("K")
names(all_Ks_C2_uniq)[2] <- paste("Cross-validation error")

plot(all_Ks_C1_uniq)
plot(all_Ks_C2_uniq)

pdf("Log-files/all_Ks_C1_uniq.pdf", height = 5, width = 5)
p1 <- ggplot(all_Ks_C1_uniq, aes(x = K, y = `Cross-validation error`)) +
    geom_point() +
    theme_classic() +
    geom_point(size = 2) +
    geom_line()
p1
dev.off()

pdf("Log-files/all_Ks_C2_uniq.pdf", height = 5, width = 5)
p2 <- ggplot(all_Ks_C2_uniq, aes(x = K, y = `Cross-validation error`)) +
    geom_point() +
    theme_classic() +
    geom_point(size = 2) +
    geom_line()
p2
dev.off()

?ggarrange

p3 <- ggarrange(p1, p2, ncol = 2, nrow = 1, labels = "AUTO")
p3

pdf("../Log-files/C1_C2_plot.pdf", height = 5, width = 10)
p3
dev.off()
