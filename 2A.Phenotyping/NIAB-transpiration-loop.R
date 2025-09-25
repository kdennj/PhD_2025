library("readxl")
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)
library(plyr)
library(tidyverse)
library("ggstatsplot")
library(rstatix)
library(ggpubr)
library(xts)
library(ggplot2)
library(MASS)
library(purrr)
library(Hmisc)
library(chron)
library(lubridate)
library(PerformanceAnalytics)
library("dtwclust") # cluster time series with dynamic time warping
library(ggdendro) # dendrograms
library(gplots) # heatmap
library(tseries) # bootstrap
library("RColorBrewer")
library(TSclust) # cluster time series

setwd("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Bean project/23.03-Field-trial-NIAB/Licor/")

#### Load all the merged data sets -----------------------------------------------

hclus_bbch <- read.xlsx("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Bean project/23.03-Field-trial-NIAB/Plots/BBCH/Clustering/DTW/DTWARP_hclus_k8_all.xlsx")
#hclus_bbch <- separate_wider_delim(hclus_bbch, cols = Accession, delim = "_", names = c("Accession", "Treat"))
hclus_bbch <- hclus_bbch[, c(1:2)]

Drought_all <- read_excel("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Bean project/23.03-Field-trial-NIAB/Plots/BBCH/Mean_all_BBCH.xlsx")
My_classification <- Drought_all[, c(1, 23)]

### Select phenotype to focus on  -----

all <- read.xlsx("E_all.xlsx")

all$first <- as.numeric(all$first)
all$second <- as.numeric(all$second)
all$third <- as.numeric(all$third)
all$fourth <- as.numeric(all$fourth)
all$fifth <- as.numeric(all$fifth)

#Square all columns 
#all <- all %>%  mutate(across(c(first, second, third, fourth, fifth), function(x) x^2))

#Scale all columns
#all <- all %>% mutate(across(c(first, second, third, fourth, fifth), function(x) scale(x)))

sum_stats <- all %>%
  group_by(Treat, Accession) %>%
  get_summary_stats(first, type = "mean_sd")
sum <- sum_stats[, c(1,2,5)] #remove sd/6 for clustering
colnames(sum) <- c("Accession", "Treat", "first", "first_sd")
sum_stats <- all %>%
  group_by(Treat, Accession) %>%
  get_summary_stats(second, type = "mean_sd")
sum_stats <- sum_stats[, c(1,2,5)]
colnames(sum_stats) <- c("Accession", "Treat", "second", "second_sd")
sum <- merge(sum, sum_stats)
sum_stats <- all %>%
  group_by(Treat, Accession) %>%
  get_summary_stats(third, type = "mean_sd")
sum_stats <- sum_stats[, c(1,2,5)]
colnames(sum_stats) <- c("Accession", "Treat", "third", "third_sd")
sum <- merge(sum, sum_stats)
sum_stats <- all %>%
  group_by(Treat, Accession) %>%
  get_summary_stats(fourth, type = "mean_sd")
sum_stats <- sum_stats[, c(1,2,5)]
colnames(sum_stats) <- c("Accession", "Treat", "fourth", "fourth_sd")
sum <- merge(sum, sum_stats)
sum_stats <- all %>%
  group_by(Treat, Accession) %>%
  get_summary_stats(fifth, type = "mean_sd")
sum_stats <- sum_stats[, c(1,2,5)]
colnames(sum_stats) <- c("Accession", "Treat", "fifth", "fifth_sd")
sum <- merge(sum, sum_stats)

sum$Accession <- paste(sum$Accession, sum$Treat, sep = "_")
df <- sum

df <- filter(df, Treat != 'C') #remove controls 

### Clustering ----
df <- df[, c(1, 3:7)]
df <- data.frame(df, row.names = 1)
colnames(df) <- c("27.07.23", "10.08.23", "18.08.23", "23.08.23", "30.08.23")

df_t <- t(df)
df_t <- as.data.frame(df_t)

models <- c("ACF", "AR.LPC.CEPS", "AR.PIC", "CDM" , "CID", "COR", "CORT" ,
            "INT.PER", "NCD", "PACF", "PDC" , "PER", "DTWARP")

#AR.MAH, DWT, PDC, PRED, MINDIST.SAX, SPEC.LLR, SPEC.GLK, SPEC.ISD failed with E

#models <- c("PER", "DTWARP")

models <- c("FRECHET")

k_options <- c(3,4,6,8,10)

for (model in models) {
 
  print(model)
  dist_ts <- TSclust::diss(SERIES = df_t, METHOD = model) # note the dataframe must be transposed
  hc <- stats::hclust(dist_ts, method="average") # method can be also "average" or diana (for DIvisive ANAlysis Clustering)
  
  for (k in k_options){
    # k for cluster which is 2 in our case (classic vs. wall)
    hclus <- stats::cutree(hc, k = k) %>% # hclus <- cluster::pam(dist_ts, k = 2)$clustering has a similar result
      as.data.frame(.) %>%
      dplyr::rename(.,cluster_group = .) %>%
      tibble::rownames_to_column("Accession") ### Use this for cluster groups information
    
    k_model <- paste(k, model, sep = "_")
    
    name <- paste("Clustering/clustering_methods/E_hclus_",k_model, ".xlsx", sep = "")
    print(name)
    
    #write.xlsx(hclus, name)
    
    df <- merge(My_classification, hclus)
    
    My_class <- df$My_Classification_less
    Clus <- df$cluster_group
    
    x <- cluster.evaluation(My_class, Clus)
    print(x)
  
}
  
}

