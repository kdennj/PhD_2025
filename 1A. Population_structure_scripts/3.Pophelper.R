#install.packages("usethis", verbose=TRUE)
#install.packages('devtools',dependencies=T)
library(devtools)
#install_github('royfrancis/pophelper')
library(pophelper)

#packageDescription("pophelper", fields="Version")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(tidyr)
library("RColorBrewer")
library(remotes)
#install.packages("reshape2")
library(reshape)
library(reshape2) 
library(plyr)
library(dplyr)

pcaPalette <- c("#2121D9","#9999FF","#FFFB23","#FF9326","#A945FF","#0089B2")
safe_pal_8 <- c("#000000", "#888888", "#5EC962", "#117733", "#332288", "#6699CC", "#CC5500","#AA4499")
match_PCA <- c("#5EC962", "#117733", "#332288", "#6699CC", "#CC5500","#AA4499")
#scales::show_col(match_PCA)


# install dependencies
#install.packages(c("ggplot2","gridExtra","label.switching","tidyr","remotes"),repos="https://cloud.r-project.org")
rm(list = ls())

library(readxl)
Seed_phenotypes<- read_excel("22.02 Seed phenotypes/22.02-updated-phenotypes-R.xlsx")


#Set the working directory depending on which ADMIXTURE files you want to analyse
#In this case I ran ADMIXTURE analysis on whole diversity panel and on the Andean 
setwd("STATS/22.01 Variant calling - GATK/22.11-ADMIXTURE-LD0.5/Q_files")

setwd("STATS/22.01 Variant calling - GATK/22.04-filtered_GATK/22.12-ADMX-K=2-vcfs/Q_files")

### ReadQ #################
#Import Qfiles
afiles <- list.files("All_Q_files", full.names=T)
afiles <- list.files("K6_Q_files/Q_files/", full.names = T)
afiles <- list.files("C2/All_Q_files", full.names = T)
afiles <- list.files("C1/All_Q_files/", full.names = T)
afiles

#Create file in correct order for grp label plots 
# Add names of accessions in this case
inds <- read.delim("samples.txt", header=FALSE, stringsAsFactors = F)
inds <- read.csv("C2/K2_Cluster_2_sample_names.csv", header=FALSE, stringsAsFactors = F)
inds <- read.csv("C1/K2_Cluster_1_sample_names.csv" , header=FALSE, stringsAsFactors = F)

Seed_phenotypes2 <- Seed_phenotypes[,c(1,3,2,4)]

#Check all sample names are correct
inds$V1[45] <- 'G50516I' #For all samples
inds$V1[29] <- 'G50516I' #For C2 ADMIXTURE
#inds$V1 <- x$V1

names(Seed_phenotypes2)[names(Seed_phenotypes2) == 'Sample_name'] <- 'V1'

#Combine the individual names and the seed phenotypes
All <- join(inds, Seed_phenotypes2)

#Add any extra information e.g. in this case country of origin
countries <- as.data.frame(All$Country)
names(countries)[names(countries) == 'All$Country'] <- 'Country'
countries_types <- All[, c(2,3)]
countries_Liborino <- All[, c(2,4)]

#Add tags to sample names for plots 
Seed_phenotypes2 <- Seed_phenotypes[,c(1,3,2,4,19)]
inds <- read.delim("samples_TAG2.txt",header=FALSE, stringsAsFactors = F)
inds$V1[45] <- 'G50516I'
Seed_phenotypes$Sample_tag <- sort(inds$V1)
names(Seed_phenotypes2)[names(Seed_phenotypes2) == 'Sample_tag'] <- 'V1'
All <- join(inds, Seed_phenotypes2)
countries_types <- All[, c(2,3)]
countries <- as.data.frame(All$Country)
Tags <- as.data.frame(All[, 1])
names(Tags)[names(Tags) == 'All[, 1]'] <- 'Tags'
Liborino <- as.data.frame(All[, 5])

#Load in file to change the order using pophelper
#Try for ADMIXTURE all samples, K6 order over 0.7
K6_clustering <- read_excel("K6_Q_files/K6_clustering_PCA-phenotypes over 0.7 ADMX.xlsx")
K6_clustering <- K6_clustering[, c(2,3,22,33)]
names(K6_clustering)[names(K6_clustering) == 'Sample_name'] <- 'V1'
#names(K6_clustering)[names(K6_clustering) == 'Sample_tag'] <- 'V1'
All2 <- join(All, K6_clustering, by = "V1" )
K6_Clus <- as.data.frame(All2$`Admixed_K6_C`)
names(K6_Clus)[names(K6_Clus) == 'All2$Admixed_A'] <- 'K6 Clusters'

K6_countries <- All2[, c(7,2)]
names(K6_countries)[names(K6_countries) == 'Admixed_A'] <- 'K6 Clusters'


#### Import K files 

K2files <- list.files("K2_Q_files", full.names = T)
K6files <- list.files("K6_Q_files/Q_files", full.names = T)
K4files <- list.files("K4_Q_files/Q_files", full.names = T)

K2files <- list.files("C1/K2_Q_files", full.names = T)
K4files <- list.files("C1/K4_Q_files", full.names = T)
K5files <- list.files("C2/K5_Q_files/Q_files", full.names = T)
K7files <- list.files("C2/K7_Q_files/Q_files", full.names = T)
K8files <- list.files("C1/K8_Q_files", full.names = T)

alist <- readQ(files=afiles)
alist <- readQ(files=K6files)
alist
class(alist)
attributes(alist)

K6list <- readQ(files = K6files)
?readQ

# add indlab to one run
rownames(alist[[1]]) <- inds$V1
# if all runs are equal length, add indlab to all runs
if(length(unique(sapply(alist,nrow)))==1) alist <- lapply(alist,"rownames<-",inds$V1)
# show row names of all runs and all samples
lapply(alist, rownames)[1:2]

# add indlab to one run
rownames(K6list[[1]]) <- inds$V1
# if all runs are equal length, add indlab to all runs
if(length(unique(sapply(K6list,nrow)))==1) K6list <- lapply(K6list,"rownames<-",inds$V1)
# show row names of all runs and all samples
lapply(K6list, rownames)[1:2]

#add tags to one run 
rownames(K7list[[1]]) <- Tags$Tags
# if all runs are equal length, add indlab to all runs
if(length(unique(sapply(K7list,nrow)))==1) K7list <- lapply(K7list,"rownames<-",Tags$Tags)
# show row names of all runs and all samples
lapply(K7list, rownames)[1:2]


##### TabulateQ ##############
tr1 <- tabulateQ(qlist = alist)
tr1

K2tr <- tabulateQ(qlist = K2list)
K4tr <- tabulateQ(qlist = K4list)
K6tr <- tabulateQ(qlist = K6list)

#### Summarise #########
sr1 <- summariseQ(tr1)
sr1

K6sr <- summariseQ(K6tr)
K6sr

K4sr <- summariseQ(K4tr)
K4sr

#######AlignK ########
alist1 <- alignK(alist, type = "across")
?alignK

K6list <- alignK(alist[51:60])
K4list <- alignK(alist[31:40])
K2list <- alignK(alist[11:20])
alist1 <- alignK(alist, type = "within")
alist2 <- alignK(alist1, type = "across")

write.csv(alist1, file = "K6_Q_files/K6_aligned2.txt")

?plotQ
c <- c(1,5,10,12,14)
plotQ(alist1[c], imgoutput = "join", exportpath=getwd(), sortind = "all", showindlab = T, sharedindlab = F)


K6align <- alignK(K6list)
write.csv(K7align, file = "C2/K7_Q_files/K7_aligned.txt")
plotQ(K6align, imgoutput = "join", exportpath=getwd(), sortind = "all", showindlab = T, sharedindlab = F)

##### Merge #########
?mergeQ
merge <- mergeQ(alist)
merge <- mergeQ(alist1)

mergeK6 <- mergeQ(K6list)
mergeK2 <- mergeQ(K2list)
mergeK4 <- mergeQ(K4list)

write.csv(merge, file = "C2/merged.txt")

#### Plotting with pophelper - at this point can move to ggplot script as cleaner plots

plotQ(merge, imgoutput = "join", exportpath=getwd(), sortind = "all", showindlab = T, sharedindlab = T,
      width = 20, height = 5, showlegend = TRUE, useindlab = T,  panelspacer = 0.1
      )
?plotQ

x <- c(1,3,5) #for ADMIXTURE all
x <- c(1,5)
x <- 1 #for C1
x <- c(4,6) #for C2

plotQ(sortQ(merge)[x], imgoutput = "join", exportpath=getwd(), sortind = "all", showindlab = T, sharedindlab = F,
      width = 28, height = 7, showlegend = TRUE, useindlab = T,  panelspacer = 0.3, legendkeysize = 7, 
      legendtextsize = 7,  indlabsize = 5, splabsize = 9, imgtype = "pdf", showyaxis = T)

plotQ(sortQ(merge)[x], imgoutput = "join", exportpath=getwd(), sortind = "label", showindlab = T, sharedindlab = T,
      width = 28, height = 7, showlegend = TRUE, useindlab = T,  panelspacer = 0.3, legendkeysize = 7, 
      legendtextsize = 7,  indlabsize = 5, splabsize = 9, imgtype = "pdf", showyaxis = T)


plotQ(sortQ(merge)[x], imgoutput = "join", exportpath=getwd(), sortind = "Cluster1", showindlab = T, sharedindlab = F,
      width = 20, height = 5, showlegend = TRUE, useindlab = T,  panelspacer = 0.1)

plotQ(merge, imgoutput = "join", exportpath=getwd(), sortind = "label", showindlab = F, sharedindlab = T,
      width = 25, height = 5, showlegend = F, useindlab = T, legendkeysize = 7, 
      legendtextsize = 7,  indlabsize = 5, splabsize = 9, showdiv = TRUE, grplab = K6_countries, 
      grplabangle = 50, ordergrp = TRUE , imgtype = "pdf", linesize = 0.4, showyaxis = T, 
      grplabheight = 1, grplabsize = 1, pointsize = 1.5, grplabpos = 0.5, linepos = 0.8, 
      panelspacer = 0.3, clustercol = match_PCA
)

###For K2, K4 and K6
write.csv(K7merge, file = "C2/K7_Q_files/K7_Clusters.txt")

K6 <- plotQ(mergeK6, returnplot = T, exportplot = F, sortind = "label", showindlab = F,
      width = 25, height = 5, showlegend = F, useindlab = T, legendkeysize = 7, 
      legendtextsize = 7,  indlabsize = 5, splabsize = 9, showdiv = TRUE, grplab = K6_countries, 
      grplabangle = 50, ordergrp = TRUE , clustercol = match_PCA, linesize = 0.4, showyaxis = T, 
      grplabheight = 1 , grplabsize = 1, pointsize = 1.5 , grplabpos = 0.5, linepos = 0.8, 
      panelspacer = 0.3
)

K2 <- plotQ(mergeK2, returnplot = T, exportplot = F, sortind = "label", showindlab = F,
      width = 25, height = 5, showlegend = F, useindlab = T, legendkeysize = 7, 
      legendtextsize = 7,  indlabsize = 5, splabsize = 9, showdiv = TRUE, grplab = K6_countries, 
      grplabangle = 50, ordergrp = TRUE , clustercol = match_PCA, linesize = 0.4, showyaxis = T, 
      grplabheight = 1 , grplabsize = 1, pointsize = 1.5 , grplabpos = 0.5, linepos = 0.8, 
      panelspacer = 0.3
)

pdf(file = "K2K4K6-ind.pdf")
grid.arrange(K2$plot[[1]], K4$plot[[1]], K6$plot[[1]], nrow=3)
dev.off()

plotQ(K6merge, exportpath=getwd(), sortind = "all", showindlab = T, useindlab = T, sharedindlab = F, width = 20, height = 8,
      showlegend = TRUE, panelspacer = 0.1, showdiv = TRUE
      )

plotQ(K4merge, exportpath=getwd(), sortind = "all", showindlab = T, useindlab = T, sharedindlab = F, width = 20, height = 8,
      showlegend = TRUE, panelspacer = 0.1, showdiv = TRUE, grplab = countries_types, grplabangle = 30, ordergrp = TRUE, 
)

write.csv(K8merge, file = "C1/K8_Q_files/K8_Clusters.txt")
plotQ(K8merge, exportpath=getwd(), sortind = "all", showindlab = T, useindlab = T, sharedindlab = F, width = 20, height = 8,
      showlegend = TRUE, panelspacer = 0.1, showdiv = TRUE
)

###Plot with gg plot ------------------------------------------------------------------------------------
#plot K6
K6_clustering <- read_excel("K6_Q_files/K6_clustering_PCA-phenotypes over 0.7 ADMX.xlsx")
K6_clustering <- K6_clustering[, c(2,26:31, 3, 33)]

K6_clustering$Sample_name <- factor(K6_clustering$Sample_name, levels = K6_clustering$Sample_name)

#### Re-order the files for plotting 
K6_clustering_order <- K6_clustering[order(K6_clustering$`Admixed_K6_C`),]
#Order <- K6_clustering$Sample_name

K6rg <- gather(K6_clustering_order, key, value, -Sample_name, -K6_order_over_0.7, -Admixed_K6_C)
#K6_rg <- K6rg[order(K6rg$K6_order_over_0.7),] 

b <- as_data_frame(K6rg$Sample_name)
b <- as.vector(b)

#plot pop proportion
cog<-ggplot(data=K6rg, aes(x=Sample_name, y=value, fill=key)) + #level = b
  geom_bar(stat="identity", show.legend = FALSE)
cog3 <- cog +theme(axis.text.x=element_text(angle = 55, hjust = 1, size=3)) +
  xlab ("Genotype") + ylab("Q value") + scale_y_continuous(expand = c(0,0)) 
cog4 <- cog3 + scale_fill_manual(values=match_PCA) + geom_vline(xintercept = K6rg$Admixed_K6_C)
#cog5 <- facet_wrap(~Admixed_K6_C, nrow = 3)

pdf(file = "plot_clusterK6.pdf",width=10,height=5)
cog4
dev.off()



#####Arrange C1/ C2 plots -----------------------------------------------------------------------------------
x <- 1 #for C1
x <- c(4,6) #for C2

#Cluster 1
p1 <- plotQ(sortQ(merge)[x], sortind = "all", showindlab = T,
            width = 28, height = 7, showlegend = TRUE, useindlab = T,  panelspacer = 0.3, legendkeysize = 7, 
            legendtextsize = 7,  indlabsize = 5, splabsize = 9, imgtype = "pdf", showyaxis = T, 
            returnplot=T,exportplot=F
            )


#Cluster 2
plotQ(merge, showindlab = F,exportpath=getwd(), sortind = "all", showindlab = T, 
            width = 25, height = 5, showlegend = F, useindlab = T,  panelspacer = 0.3, legendkeysize = 7, 
            legendtextsize = 7,  indlabsize = 5, splabsize = 9, imgtype = "pdf", showyaxis = T, 
            returnplot=F,exportplot=T, showdiv = TRUE
            )
?plotQ

plotQ(sortQ(merge)[x], imgoutput = "join", exportpath=getwd(), sortind = "label", showindlab = F, sharedindlab = T,
      width = 25, height = 5, showlegend = F, useindlab = T, legendkeysize = 7, 
      legendtextsize = 7,  indlabsize = 5, splabsize = 9, showdiv = TRUE, grplab = K6_countries, 
      grplabangle = 50, ordergrp = TRUE , imgtype = "pdf", linesize = 0.4, showyaxis = T, 
      grplabheight = 1, grplabsize = 1, pointsize = 1.5, grplabpos = 0.5, linepos = 0.8, 
      panelspacer = 0.3, clustercol = match_PCA
)

p3 <- grid.arrange(p1$plot[[1]], p2$plot[[1]], ncol = 1, nrow = 2, heights = c(1,2))
p3

pdf("C1K2_C2K5K7.pdf", height = 5, width = 5)
p3
dev.off()

#?plotQMultiline
#plotQMultiline(K6merge, sortind = "all", showindlab = T, exportpath=getwd())

#plotQ(readQ(list.files(pattern = "merged", full.names = T)), imgoutput = "join", sortind = "all", showindlab = T, sharedindlab = F, basesize = 10,width = 30, height = 10)
#collectClumppOutput(filetype="merged")
