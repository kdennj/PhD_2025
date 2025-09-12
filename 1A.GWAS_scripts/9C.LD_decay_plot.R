rm(list = ls())


#############Plot from PopLDdecay -------------------------
library(scales)
library(dplyr)
library(stringr)
library(ggplot2)

setwd("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Bean project/STATS/22.01 Variant calling - GATK/popLDdecay")

ld <- read.delim("LDdecay.stat", stringsAsFactors = FALSE,header=TRUE, sep = "\t") 
#ld <- read.delim("admixture/LDdecay.stat", stringsAsFactors = FALSE,header=TRUE, sep = "\t") 
#ld <- read.delim("LDdecay_Chr01.stat", stringsAsFactors = FALSE,header=TRUE, sep = "\t") 
ld <- read.delim("LDdecay_M.stat", stringsAsFactors = FALSE,header=TRUE, sep = "\t") 


meanld <- mean(ld$Mean_r.2)

maxld <- max(ld$Mean_r.2) #using max LD value from adjusted data
halfdecay = maxld*0.5
halfdecay <- 0.24
halfdecaydist <- ld$X.Dist[which.min(abs(ld$Mean_r.2-halfdecay))]

jpeg("PopLDdecay.jpg", width = 400, height = 400)
p <- ggplot(ld, aes(X.Dist, Mean_r.2)) + geom_point(size = 1) +
  theme_classic() +
  xlab("Distance (bp)") + ylab("Mean r²") +
  geom_smooth() +
  geom_hline(yintercept = 0.1, colour = "red") +
  geom_vline(xintercept = halfdecaydist, colour = "darkgreen") +
  scale_y_continuous(limits = c(0, 1), labels = label_number()) 
p
dev.off()

