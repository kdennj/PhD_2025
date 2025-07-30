library(pophelper)
packageDescription("pophelper", fields="Version")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(tidyr)
library(remotes)
library(plyr)
library(dplyr)

rm(list = ls())

match_PCA <- c("#5EC962", "#117733", "#332288", "#6699CC", "#CC5500","#AA4499")

library(readxl)
Seed_phenotypes<- read_excel("Bean project/Common bean samples/22.02 Seed phenotypes/22.02-updated-phenotypes-R.xlsx")

setwd("/Bean project/STATS/22.01 Variant calling - GATK/23.01-STRUCTURE/")

#Create file in correct order for grp label plots 
#x <- c(1:144)
inds <- read.delim("samples.txt", header=FALSE, stringsAsFactors = F)

Seed_phenotypes2 <- Seed_phenotypes[,c(1,3,2)]
inds$V1[45] <- 'G50516I' #For all samples
#inds$V1 <- x$V1
names(Seed_phenotypes2)[names(Seed_phenotypes2) == 'Sample_name'] <- 'V1'
All <- join(inds, Seed_phenotypes2)
countries <- as.data.frame(All$Country)
names(countries)[names(countries) == 'All$Country'] <- 'Country'
countries_types <- All[, c(2,3)]

#Add tags to sample names for plots 
Seed_phenotypes2 <- Seed_phenotypes[,c(1,3,2,4,19)]
inds <- read.delim("samples_TAG2.txt",header=FALSE, stringsAsFactors = F)
inds$V1[45] <- 'G50516I'
Seed_phenotypes$Sample_tag <- sort(inds$V1)
names(Seed_phenotypes2)[names(Seed_phenotypes2) == 'Sample_tag'] <- 'V1'
All <- join(inds, Seed_phenotypes2)
countries_types <- All[, c(2,3)]
Tags <- as.data.frame(All[, 1])
names(Tags)[names(Tags) == 'All[, 1]'] <- 'Tags'
countries <- as.data.frame(All$Country)
Liborino <- as.data.frame(All[, 5])

#Load in file to change the order using pophelper
#Try for ADMIXTURE all samples, K6 order over 0.7
K6_clustering <- read_excel("/Bean project/STATS/22.01 Variant calling - GATK/22.11-ADMIXTURE-LD0.5/Q_files/K6_Q_files/K6_clustering_PCA-phenotypes over 0.7 ADMX.xlsx")
K6_clustering <- K6_clustering[, c(2,3,22,33)]
names(K6_clustering)[names(K6_clustering) == 'Sample_name'] <- 'V1'
#names(K6_clustering)[names(K6_clustering) == 'Sample_tag'] <- 'V1'
All2 <- join(All, K6_clustering, by = "V1" )
K6_Clus <- as.data.frame(All2$Admixed_A)
names(K6_Clus)[names(K6_Clus) == 'All2$Admixed_A'] <- 'K6 Clusters'

K6_countries <- All2[, c(6,2)]
names(K6_countries)[names(K6_countries) == 'Admixed_A'] <- 'K6 Clusters'

### ReadQ #################
afiles <- list.files("K6", full.names=T)
afiles <- list.files("Results", full.names=T)
afiles

alist <- readQ(files=afiles)
alist
class(alist)
attributes(alist)

# add indlab to one run
rownames(alist[[1]]) <- inds$V1
# if all runs are equal length, add indlab to all runs
if(length(unique(sapply(alist,nrow)))==1) alist <- lapply(alist,"rownames<-",inds$V1)
# show row names of all runs and all samples
lapply(alist, rownames)[1:2]


##### TabulateQ ##############
tr1 <- tabulateQ(qlist = alist)
tr1

write.csv(tr1,"K6_tablulateQ.csv")

#### Summarise #########
sr1 <- summariseQ(tr1)
sr1

write.csv(Sr1,"K6_SummariseQ.csv")

####Evanno Method for STRUCTURE#####
evannoMethodStructure(data=sr1,exportplot=T,writetable=T,na.rm=T,exportpath=getwd())
?evannoMethodStructure

###Align########
alist1 <- alignK(alist, type = "across")
?alignK

alist1 <- alignK(alist, type = "within")
alist2 <- alignK(alist1, type = "across")

write.csv(alist1, file = "aligned.txt")
plotQ(alist1, imgoutput = "join", exportpath=getwd(), sortind = "all", showindlab = T, sharedindlab = F)

##### Merge #########
?mergeQ
merge <- mergeQ(alist2)

write.csv(merge, file = "merged.txt")

plotQ(merge, imgoutput = "join", exportpath=getwd(), sortind = "all", showindlab = T, sharedindlab = F,
      width = 20, height = 5, showlegend = TRUE, useindlab = T,  panelspacer = 0.1
)
?plotQ

x <- c(1,5)

plotQ(sortQ(merge)[x], imgoutput = "join", exportpath=getwd(), sortind = "all", showindlab = T, sharedindlab = F,
      width = 28, height = 10, showlegend = TRUE, useindlab = T,  panelspacer = 0.2, legendkeysize = 7, 
      legendtextsize = 7,  indlabsize = 5, splabsize = 9, imgtype = "pdf", showyaxis = T)

plotQ(sortQ(merge)[x], imgoutput = "join", exportpath=getwd(), sortind = "label", showindlab = T, sharedindlab = T,
      width = 28, height = 7, showlegend = TRUE, useindlab = T,  panelspacer = 0.3, legendkeysize = 7, 
      legendtextsize = 7,  indlabsize = 5, splabsize = 9, imgtype = "pdf", showyaxis = T)


plotQ(merge[c(5)], exportpath=getwd(), sortind = "all", showindlab = T, sharedindlab = F,
      width = 20, height = 5, showlegend = TRUE, useindlab = T,  panelspacer = 0.1
)


plotQ(merge, exportpath=getwd(), sortind = "all", showindlab = T, useindlab = T, sharedindlab = F, width = 20, height = 8,
      showlegend = TRUE, panelspacer = 0.1, showdiv = TRUE
)


plotQ(sortQ(merge)[x], imgoutput = "sep", exportpath=getwd(), showindlab = F, sharedindlab = F,
      width = 25, height = 10, showlegend = TRUE, useindlab = T, legendkeysize = 7, 
      legendtextsize = 7,  indlabsize = 5, splabsize = 9, showdiv = TRUE, grplab = countries, 
      grplabangle = 50, ordergrp = TRUE , imgtype = "pdf", linesize = 0.4, showyaxis = T, 
      grplabheight = 1, grplabsize = 1.8, pointsize = 1.5, grplabpos = 0.5, linepos = 0.8, 
      panelspacer = 0.3
)

plotQ(sortQ(merge)[x], imgoutput = "join", exportpath=getwd(), sortind = "label", showindlab = F, sharedindlab = T,
      width = 25, height = 5, showlegend = T, useindlab = T, legendkeysize = 7, 
      legendtextsize = 7,  indlabsize = 5, splabsize = 9, showdiv = TRUE, grplab = K6_countries, 
      grplabangle = 50, ordergrp = TRUE , imgtype = "pdf", linesize = 0.4, showyaxis = T, 
      grplabheight = 1, grplabsize = 1, pointsize = 1.5, grplabpos = 0.5, linepos = 0.8, 
      panelspacer = 0.3, clustercol = match_PCA
)

