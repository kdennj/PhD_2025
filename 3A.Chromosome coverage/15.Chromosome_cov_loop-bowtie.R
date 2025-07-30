library("RColorBrewer")
library("ggplot2")

library(dplyr)
library(tidyr)
library(stringr)
library("gridExtra")
library("cowplot")
library(patchwork)
library("ggpubr")


#Following what Janet used - median
setwd("/STATS/22.03 Chromosome coverage/bowtie2/23.02-bedgraphs-mean/")

rm(list = ls())

# normalise by mean
# calculate mean coverage for A and B genome - store in list and output as dataframe
# calculate mean coverage for A and B genome per chromsome - output as csv
# do AB plots for each sample

# read in list of samples 
pop <- read.table("../../samples_uniq_all.txt") #list of sample names

slist <- pop$V1
str(slist)

 # initialise empty lists
mean_list <- c()
normA_list <- c()
normM_list <- c()
Amean_list <- c()
Mmean_list <- c()

# read in coverage file for each sample
# output mean coverage over each chromosoem as csv
# plot normalised coverage by A genome
# plot normalised coverage by B genome

for (sam in slist) {
  name <- paste(sam,"_unique.bam.bedgraph_median_cov.txt",sep = "")
  print (name) 
  #all <- read.delim ("S0001_G50516H_unique.bam_rep.bedgraph_mean_cov_rep.txt", header = FALSE)
  all<- read.delim (name, header = FALSE)
  colnames(all) <- c("chrom","start","end","cov")
  head(all)
  
  # to set window with no coverage of 0, these are input as "." 
  all$cov <- as.numeric(all$cov) 
  all$cov[is.na(all$cov)] <- 0
  all$cov <- as.numeric(all$cov) 
  print (sam)
  #colnames(all) <- c("chrom","start","end","cov")
  cov_chr <- all %>% group_by(chrom) %>%
    summarize(cov = mean(cov))
  ccname <- paste(sam,"_chr_cov.csv",sep = "")
  write.csv(cov_chr,ccname)
  ma <- mean(all$cov)
  print(ma)
  mean_list <- append(mean_list,c(ma))
  all$covNormByMean <- all$cov/mean(all$cov)  ###relative coverage
  check <- mean(all$covNormByMean)
  print(check) #Should =1
  # get A genome chromosomes
  Achr01<- all %>% filter(grepl('^ChrA', chrom))
  maA <-  mean(Achr01$cov)
  normA_list <- append(normA_list,c(maA))
  maAn <- mean(Achr01$covNormByMean)
  Amean_list <- append(Amean_list,c(maAn))
  # get M genome chromosomes
  Mchr01<- all %>% filter(grepl('^ChrM', chrom))
  maM <- mean(Mchr01$cov)
  normM_list <- append(normM_list,c(maM))
  maMn <- mean(Mchr01$covNormByMean)
  Mmean_list <- append(Mmean_list,c(maMn))
  
  ##Define function 
  all.neg <- function(x) -1*abs(x)
  
  ## make B negative and remove B from chromosome -
  
  #Mchr01neg <- cbind(Mchr01,"covbynorm_neg" = all.neg(Mchr01$covNormByMean))
  Mchr01neg <- cbind(Mchr01,"cov_neg" = all.neg(Mchr01$cov))
  
  #Achr01neg <- cbind(Achr01,"covbynorm_neg" = all.neg(Achr01$covNormByMean))
  Achr01neg <- cbind(Achr01,"cov_neg" = all.neg(Achr01$cov))
  
  Mnes <- Mchr01neg %>%
    mutate(chrom = str_replace(chrom,"^ChrM","Chr"))
  
  Mpos <- Mchr01 %>%
    mutate(chrom = str_replace(chrom,"^ChrM","Chr"))
  
  Apos <- Achr01 %>% 
    mutate(chrom = str_replace(chrom, "^ChrA", "Chr"))
  
  Anes <- Achr01neg %>% 
    mutate(chrom = str_replace(chrom, "^ChrA", "Chr"))
  
  #print (head(Mnes))
  #print (head(Mchr01))
  
  #######################
  ###replot with syntheny
  #######################
  #read minimap alignments: 
  minimap <- read.delim("../../M_100kbwindows-over-A.longest_alignmentotal.paf", header=F)
  minimap <- minimap[,c(1,6,8,9)]
  colnames(minimap) <- c("uniqID","inAA_chr","inAA_start","inAA_stop")
  minimap <- minimap[- grep("scaffoldA", minimap$inAA_chr),]
  #create id field in Mchr01neg
  Mchr01neg$uniqID <- paste0(Mchr01neg$chrom,":",Mchr01neg$start,"-",Mchr01neg$end)
  
  #join left both dataframes
  Mchr01negJOIN <- merge(x=Mchr01neg, y=minimap, by="uniqID", all.x = FALSE) #false; so as to ignore Ms without A homologous
  
  Achr01$facet <- Achr01$chrom
  Mchr01negJOIN$facet <- Mchr01negJOIN$inAA_chr
  
  aname <- paste(sam,"cov_byA.pdf",sep = "")
  #aname <- 'file_cov_byA'
  
  print (aname)
  print (head(Achr01))
  
  
  p <- ggplot() + geom_bar(data=Achr01,aes(x=start,y=cov), stat = "identity",color="dodgerblue")
  pg <- p + geom_bar(data=Mchr01negJOIN,aes(x=inAA_start,y=cov_neg), stat = "identity",color="firebrick3")
  pgx <- pg + coord_cartesian(ylim = c(-15,15)) + theme(axis.text = element_text(size = 5))  
  pgxl <- pgx + geom_hline(yintercept=0, linetype="solid", color = "gray40")
  pgxk <- pgxl + geom_hline(yintercept=maA, linetype="dashed", color = "gray30", alpha=0.3)
  pgxm <- pgxk + geom_hline(yintercept=-maM, linetype="dashed", color = "gray30", alpha=0.3)
  pgsw <- pgxm + facet_wrap(~facet, ncol=1, strip.position = "right")
  pgs <- pgsw + theme_classic()
  pgs
  
  pdf(file = paste("plots_by_cov/", aname,sep = ""), width =8, height = 10)
  plot(pgs)
  dev.off()
  
  
  ###SAME WITH M AS REFERENCE (ALSO SWAPPING TOP AND BOTTOM BUT NOT COLOURS):
  minimapM <- read.delim("../../A_100kbwindows-over-M.longest_alignmentotal.paf", header=F)
  minimapM <- minimapM[,c(1,6,8,9)]
  colnames(minimapM) <- c("uniqID","inM_chr","inM_start","inM_stop")
  minimapM <- minimapM[- grep("scaffoldM", minimapM$inM_chr),]
  #create id field in Bchr01neg
  Achr01neg$uniqID <- paste0(Achr01neg$chrom,":",Achr01neg$start,"-",Achr01neg$end)
  
  #join left both dataframes
  Achr01negJOIN <- merge(x=Achr01neg, y=minimapM, by="uniqID", all.x = FALSE) #false; so as to ignore As without M homologous
  
  Achr01negJOIN$facet <- Achr01negJOIN$inM_chr
  Mchr01$facet <- Mchr01$chrom
  
  bname <- paste(sam,"cov_byM.pdf",sep = "")
  #bname <- 'file_cov_byM'
  
  p <- ggplot() + geom_bar(data=Mchr01,aes(x=start,y=cov), stat = "identity",color="firebrick3")
  pg <- p + geom_bar(data=Achr01negJOIN,aes(x=inM_start,y=cov_neg), stat = "identity",color="dodgerblue")
  pgx <- pg + coord_cartesian(ylim = c(-15,15)) + theme(axis.text = element_text(size = 5))  
  pgxl <- pgx + geom_hline(yintercept=0, linetype="solid", color = "gray40")
  pgxk <- pgxl + geom_hline(yintercept=maM, linetype="dashed", color = "gray30", alpha=0.3)
  pgxm <- pgxk + geom_hline(yintercept=-maA, linetype="dashed", color = "gray30", alpha=0.3)
  pgsw <- pgxm + facet_wrap(~facet,ncol=1,strip.position = "right")
  pgsM <- pgsw + theme_classic()
  pgsM
  
  pdf(file = paste("plots_by_cov/", bname,sep = ""),width =8, height = 10)
  plot(pgsM)
  dev.off()
  
  ###Combine both plots
  p5 = ggarrange(pgs, pgsM, ncol = 2, nrow = 1, legend = "none")
  p5
  pdf(file = paste("plots_grouped_together/", aname,sep = ""), width =12, height = 12)
  
  plot(p5)
  dev.off()
  
}



# make dataframe from lists and output as csv
new_df <- as.data.frame(slist)
results <- cbind(new_df,mean_list,normA_list,Amean_list,normM_list,Mmean_list)
write.csv(results,"AM_meanAM_results.csv")

# check lists are all the same length
length(slist)
length(mean_list)
length(normA_list)
length(normM_list)
length(Amean_list)
length(Mmean_list)

#----------------------------------
####Cov by mean - normalised
p <- ggplot() + geom_bar(data=Achr01,aes(x=start,y=covbymean), stat = "identity",color="dodgerblue")
pg <- p + geom_bar(data=Mchr01negJOIN,aes(x=inAA_start,y=covbymean_neg), stat = "identity",color="firebrick3")
pgx <- pg + coord_cartesian(ylim = c(-5,5)) + theme(axis.text = element_text(size = 5))  
pgxl <- pgx + geom_hline(yintercept=0, linetype="solid", color = "gray40")
pgxk <- pgxl + geom_hline(yintercept=1.04, linetype="dashed", color = "gray30", alpha=0.3)
pgxm <- pgxk + geom_hline(yintercept=-0.92, linetype="dashed", color = "gray30", alpha=0.3)
pgsw <- pgxm + facet_wrap(~facet, ncol=1, strip.position = "right")
pgs <- pgsw + theme_classic()
pgs

pdf(file = aname,width =8, height = 12)
plot(pgs)
dev.off()
  
########
#Following what Janet used - median
#I prefer median as includes all data 
setwd("/STATS/22.03 Chromosome coverage/bowtie2/23.02-bedgraphs-mean/")

rm(list = ls())

# normalise by mean
# calculate mean coverage for A and B genome - store in list and output as dataframe
# calculate mean coverage for A and B genome per chromsome - output as csv
# do AB plots for each sample

# read in list of samples 
pop <- read.table("../../bowtie_practise.txt") #list of sample names

slist <- pop$V1
str(slist)


# initialise empty lists
mean_list <- c()
normA_list <- c()
normM_list <- c()
Amean_list <- c()
Mmean_list <- c()

for (sam in slist) {
  name <- paste(sam,"_unique.bam.bedgraph_mean_cov.txt",sep = "")
  print (name) 
  X <- sam
  #all <- read.delim ("S0001_G50516H_pp.bam_rep.bedgraph_median_cov.txt", header = FALSE)
  all<- read.delim (name, header = FALSE)
  colnames(all) <- c("chrom","start","end","cov")
  head(all)
  
  # to set window with no coverage of 0, these are input as "." 
  all$cov <- as.numeric(all$cov) 
  all$cov[is.na(all$cov)] <- 0
  
  print (sam)
  #colnames(all) <- c("chrom","start","end","cov")
  
  # cov_chr <- all %>% group_by(chrom) %>%
  #   summarise(cov = mean(cov))
  
  cov_chr <- aggregate(all[,4], list(all$chrom), mean)
  
  ccname <- paste(sam,"_chr_cov.csv",sep = "")
  write.csv(cov_chr,ccname)
  
  ma <- mean(all$cov)
  print(ma)
  mean_list <- append(mean_list,c(ma))
  
  all$covNormByMean <- all$cov/mean(all$cov)  ###relative coverage
  check <- mean(all$covNormByMean)
  print(check) #Should =1
  
  # get A genome chromosomes
  Achr01<- all %>% filter(grepl('^ChrA', chrom))
  maA <-  mean(Achr01$cov)
  normA_list <- append(normA_list,c(maA))
  maAn <- mean(Achr01$covNormByMean)
  Amean_list <- append(Amean_list,c(maAn))
  
  # get M genome chromosomes
  Mchr01<- all %>% filter(grepl('^ChrM', chrom))
  maM <- mean(Mchr01$cov)
  normM_list <- append(normM_list,c(maM))
  maMn <- mean(Mchr01$covNormByMean)
  Mmean_list <- append(Mmean_list,c(maMn))
  
  ##Define function 
  all.neg <- function(x) -1*abs(x)
  
  ## make B negative and remove B from chromosome -
  
  Mchr01neg <- cbind(Mchr01,"covbynorm_neg" = all.neg(Mchr01$covNormByMean))
  #Mchr01neg <- cbind(Mchr01,"cov_neg" = all.neg(Mchr01$cov))
  
  Achr01neg <- cbind(Achr01,"covbynorm_neg" = all.neg(Achr01$covNormByMean))
  #Achr01neg <- cbind(Achr01,"cov_neg" = all.neg(Achr01$cov))
  
  Mnes <- Mchr01neg %>%
    mutate(chrom = str_replace(chrom,"^ChrM","Chr"))
  
  Mpos <- Mchr01 %>%
    mutate(chrom = str_replace(chrom,"^ChrM","Chr"))
  
  Apos <- Achr01 %>% 
    mutate(chrom = str_replace(chrom, "^ChrA", "Chr"))
  
  Anes <- Achr01neg %>% 
    mutate(chrom = str_replace(chrom, "^ChrA", "Chr"))
  
  #print (head(Mnes))
  #print (head(Mchr01))
  
  #######################
  ###replot with syntheny
  #######################
  #read minimap alignments: 
  minimap <- read.delim("../../M_100kbwindows-over-A.longest_alignmentotal.paf", header=F)
  minimap <- minimap[,c(1,6,8,9)]
  colnames(minimap) <- c("uniqID","inAA_chr","inAA_start","inAA_stop")
  minimap <- minimap[- grep("scaffoldA", minimap$inAA_chr),]
  #create id field in Mchr01neg
  Mchr01neg$uniqID <- paste0(Mchr01neg$chrom,":",Mchr01neg$start,"-",Mchr01neg$end)
  
  #join left both dataframes
  Mchr01negJOIN <- merge(x=Mchr01neg, y=minimap, by="uniqID", all.x = FALSE) #false; so as to ignore Ms without A homologous
  
  Achr01$facet <- Achr01$chrom
  Mchr01negJOIN$facet <- Mchr01negJOIN$inAA_chr
  
  aname <- paste(sam,"cov_byA.pdf",sep = "")
  #aname <- 'file_cov_byA'
  
  print (aname)
  print (head(Achr01))
  
  ####Cov by mean - normalised
  p <- ggplot() + geom_bar(data=Achr01,aes(x=start,y=covNormByMean), stat = "identity",color="dodgerblue")
  pg <- p + geom_bar(data=Mchr01negJOIN,aes(x=inAA_start,y=covbynorm_neg), stat = "identity",color="firebrick3")
  pgx <- pg + coord_cartesian(ylim = c(-2.5,2.5)) + theme(axis.text = element_text(size = 5))  
  pgxl <- pgx + geom_hline(yintercept=0, linetype="solid", color = "gray40")
  pgxk <- pgxl + geom_hline(yintercept=maAn, linetype="dashed", color = "gray30", alpha=0.3)
  pgxm <- pgxk + geom_hline(yintercept=-maMn, linetype="dashed", color = "gray30", alpha=0.3)
  pgsw <- pgxm + facet_wrap(~facet, ncol=1, strip.position = "right")
  pgs <- pgsw + theme_classic()
  pgs
  
  pdf(file = paste("plots_by_norm-genome-mean/", aname,sep = ""), width =8, height = 12)
  plot(pgs)
  dev.off()
  
  write.csv(Achr01, file = paste(X, "_Achr.csv", sep = ""))
  #write.csv(Mchr01negJOIN, file = paste(X, "_MchrnegJOIN.csv", sep = ""))
  
  ###SAME WITH M AS REFERENCE (ALSO SWAPPING TOP AND BOTTOM BUT NOT COLOURS):
  minimapM <- read.delim("../../A_100kbwindows-over-M.longest_alignmentotal.paf", header=F)
  minimapM <- minimapM[,c(1,6,8,9)]
  colnames(minimapM) <- c("uniqID","inM_chr","inM_start","inM_stop")
  minimapM <- minimapM[- grep("scaffoldM", minimapM$inM_chr),]
  #create id field in Bchr01neg
  Achr01neg$uniqID <- paste0(Achr01neg$chrom,":",Achr01neg$start,"-",Achr01neg$end)
  
  #join left both dataframes
  Achr01negJOIN <- merge(x=Achr01neg, y=minimapM, by="uniqID", all.x = FALSE) #false; so as to ignore As without M homologous
  
  Achr01negJOIN$facet <- Achr01negJOIN$inM_chr
  Mchr01$facet <- Mchr01$chrom
  
  bname <- paste(sam,"cov_byM.pdf",sep = "")
  #bname <- 'file_cov_byM'
  
  ####Cov by mean - normalised
  p <- ggplot() + geom_bar(data=Mchr01,aes(x=start,y=covNormByMean), stat = "identity",color="firebrick3")
  pg <- p + geom_bar(data=Achr01negJOIN,aes(x=inM_start,y=covbynorm_neg), stat = "identity",color="dodgerblue")
  pgx <- pg + coord_cartesian(ylim = c(-2.5,2.5)) + theme(axis.text = element_text(size = 5))  
  pgxl <- pgx + geom_hline(yintercept=0, linetype="solid", color = "gray40")
  pgxk <- pgxl + geom_hline(yintercept=maMn, linetype="dashed", color = "gray30", alpha=0.3)
  pgxm <- pgxk + geom_hline(yintercept=-maAn, linetype="dashed", color = "gray30", alpha=0.3)
  pgsw <- pgxm + facet_wrap(~facet, ncol=1, strip.position = "right")
  pgsM <- pgsw + theme_classic()
  pgsM
  
  pdf(file = paste("plots_by_norm-genome-mean/", bname,sep = ""), width =8, height = 12)
  plot(pgsM)
  dev.off()
  
  write.csv(Mchr01, file = paste(X, "_Mchr.csv", sep = ""))
  #write.csv(Achr01negJOIN, file = paste(X, "_AchrnegJOIN.csv", sep = ""))
  
  ###Combine both plots
  p5 = ggarrange(pgs, pgsM, ncol = 2, nrow = 1, legend = "none")
  p5
  pdf(file = paste("plots_by_norm_grouped_together-genome-mean/", aname,sep = ""), width =12, height = 12)
  
  plot(p5)
  dev.off()
  
}

########
#Testing normalising mean by each chromosome not whole genome ############
setwd("/STATS/22.03 Chromosome coverage/bowtie2/23.02-bedgraphs-mean/")

rm(list = ls())

# normalise by mean
# calculate mean coverage for A and B genome - store in list and output as dataframe
# calculate mean coverage for A and B genome per chromsome - output as csv
# do AB plots for each sample

# read in list of samples 
pop <- read.table("../../samples_uniq_all.txt") #list of sample names

slist <- pop$V1
str(slist)


# initialise empty lists
mean_list <- c()
normA_list <- c()
normM_list <- c()
Amean_list <- c()
Mmean_list <- c()

for (sam in slist) {
  name <- paste(sam,"_unique.bam.bedgraph_mean_cov.txt",sep = "")
  print (name) 
  X <- sam
  #all <- read.delim ("S0001_G50516H_pp.bam_rep.bedgraph_median_cov.txt", header = FALSE)
  all<- read.delim (name, header = FALSE)
  colnames(all) <- c("chrom","start","end","cov")
  head(all)
  
  # to set window with no coverage of 0, these are input as "." 
  all$cov <- as.numeric(all$cov) 
  all$cov[is.na(all$cov)] <- 0
  
  print (sam)
  #colnames(all) <- c("chrom","start","end","cov")
  
  cov_chr <- aggregate(all[,4], list(all$chrom), mean)
  colnames(cov_chr) <- c("chrom", "mean_cov")
  
  all <- merge(all, cov_chr, by = "chrom")
  
  ccname <- paste(sam,"_chr_cov.csv",sep = "")
  write.csv(cov_chr,ccname)
  
  ma <- mean(all$cov)
  print(ma)
  mean_list <- append(mean_list,c(ma))
  
  all$covNormByMean <- all$cov/all$mean_cov ###relative coverage
  check <- mean(all$covNormByMean)
  print(check) #Should =1
  
  # get A genome chromosomes
  Achr01<- all %>% filter(grepl('^ChrA', chrom))
  maA <-  mean(Achr01$cov)
  normA_list <- append(normA_list,c(maA))
  maAn <- mean(Achr01$covNormByMean)
  Amean_list <- append(Amean_list,c(maAn))
  
  # get M genome chromosomes
  Mchr01<- all %>% filter(grepl('^ChrM', chrom))
  maM <- mean(Mchr01$cov)
  normM_list <- append(normM_list,c(maM))
  maMn <- mean(Mchr01$covNormByMean)
  Mmean_list <- append(Mmean_list,c(maMn))
  
  ##Define function 
  all.neg <- function(x) -1*abs(x)
  
  ## make B negative and remove B from chromosome -
  
  Mchr01neg <- cbind(Mchr01,"covbynorm_neg" = all.neg(Mchr01$covNormByMean))
  #Mchr01neg <- cbind(Mchr01,"cov_neg" = all.neg(Mchr01$cov))
  
  Achr01neg <- cbind(Achr01,"covbynorm_neg" = all.neg(Achr01$covNormByMean))
  #Achr01neg <- cbind(Achr01,"cov_neg" = all.neg(Achr01$cov))
  
  Mnes <- Mchr01neg %>%
    mutate(chrom = str_replace(chrom,"^ChrM","Chr"))
  
  Mpos <- Mchr01 %>%
    mutate(chrom = str_replace(chrom,"^ChrM","Chr"))
  
  Apos <- Achr01 %>% 
    mutate(chrom = str_replace(chrom, "^ChrA", "Chr"))
  
  Anes <- Achr01neg %>% 
    mutate(chrom = str_replace(chrom, "^ChrA", "Chr"))
  
  #print (head(Mnes))
  #print (head(Mchr01))
  
  #######################
  ###replot with syntheny
  #######################
  #read minimap alignments: 
  minimap <- read.delim("../../M_100kbwindows-over-A.longest_alignmentotal.paf", header=F)
  minimap <- minimap[,c(1,6,8,9)]
  colnames(minimap) <- c("uniqID","inAA_chr","inAA_start","inAA_stop")
  minimap <- minimap[- grep("scaffoldA", minimap$inAA_chr),]
  #create id field in Mchr01neg
  Mchr01neg$uniqID <- paste0(Mchr01neg$chrom,":",Mchr01neg$start,"-",Mchr01neg$end)
  
  #join left both dataframes
  Mchr01negJOIN <- merge(x=Mchr01neg, y=minimap, by="uniqID", all.x = FALSE) #false; so as to ignore Ms without A homologous
  
  Achr01$facet <- Achr01$chrom
  Mchr01negJOIN$facet <- Mchr01negJOIN$inAA_chr
  
  aname <- paste(sam,"cov_byA.pdf",sep = "")
  #aname <- 'file_cov_byA'
  
  print (aname)
  print (head(Achr01))
  
  ####Cov by mean - normalised
  p <- ggplot() + geom_bar(data=Achr01,aes(x=start,y=covNormByMean), stat = "identity",color="dodgerblue")
  pg <- p + geom_bar(data=Mchr01negJOIN,aes(x=inAA_start,y=covbynorm_neg), stat = "identity",color="firebrick3")
  pgx <- pg + coord_cartesian(ylim = c(-2.5,2.5)) + theme(axis.text = element_text(size = 5))  
  pgxl <- pgx + geom_hline(yintercept=0, linetype="solid", color = "gray40")
  pgxk <- pgxl + geom_hline(yintercept=maAn, linetype="dashed", color = "gray30", alpha=0.3)
  pgxm <- pgxk + geom_hline(yintercept=-maMn, linetype="dashed", color = "gray30", alpha=0.3)
  pgsw <- pgxm + facet_wrap(~facet, ncol=1, strip.position = "right")
  pgs <- pgsw + theme_classic()
  pgs
  
  pdf(file = paste("plots_by_norm_mean-by-chrom/", aname,sep = ""), width =8, height = 12)
  plot(pgs)
  dev.off()
  
  write.csv(Achr01, file = paste(X, "_Achr.csv", sep = ""))
  #write.csv(Mchr01negJOIN, file = paste(X, "_MchrnegJOIN.csv", sep = ""))
  
  ###SAME WITH M AS REFERENCE (ALSO SWAPPING TOP AND BOTTOM BUT NOT COLOURS):
  minimapM <- read.delim("../../A_100kbwindows-over-M.longest_alignmentotal.paf", header=F)
  minimapM <- minimapM[,c(1,6,8,9)]
  colnames(minimapM) <- c("uniqID","inM_chr","inM_start","inM_stop")
  minimapM <- minimapM[- grep("scaffoldM", minimapM$inM_chr),]
  #create id field in Bchr01neg
  Achr01neg$uniqID <- paste0(Achr01neg$chrom,":",Achr01neg$start,"-",Achr01neg$end)
  
  #join left both dataframes
  Achr01negJOIN <- merge(x=Achr01neg, y=minimapM, by="uniqID", all.x = FALSE) #false; so as to ignore As without M homologous
  
  Achr01negJOIN$facet <- Achr01negJOIN$inM_chr
  Mchr01$facet <- Mchr01$chrom
  
  bname <- paste(sam,"cov_byM.pdf",sep = "")
  #bname <- 'file_cov_byM'
  
  ####Cov by mean - normalised
  p <- ggplot() + geom_bar(data=Mchr01,aes(x=start,y=covNormByMean), stat = "identity",color="firebrick3")
  pg <- p + geom_bar(data=Achr01negJOIN,aes(x=inM_start,y=covbynorm_neg), stat = "identity",color="dodgerblue")
  pgx <- pg + coord_cartesian(ylim = c(-2.5,2.5)) + theme(axis.text = element_text(size = 5))  
  pgxl <- pgx + geom_hline(yintercept=0, linetype="solid", color = "gray40")
  pgxk <- pgxl + geom_hline(yintercept=maMn, linetype="dashed", color = "gray30", alpha=0.3)
  pgxm <- pgxk + geom_hline(yintercept=-maAn, linetype="dashed", color = "gray30", alpha=0.3)
  pgsw <- pgxm + facet_wrap(~facet, ncol=1, strip.position = "right")
  pgsM <- pgsw + theme_classic()
  pgsM
  
  pdf(file = paste("plots_by_norm_mean-by-chrom/", bname,sep = ""), width =8, height = 12)
  plot(pgsM)
  dev.off()
  
  write.csv(Mchr01, file = paste(X, "_Mchr.csv", sep = ""))
  #write.csv(Achr01negJOIN, file = paste(X, "_AchrnegJOIN.csv", sep = ""))
  
  ###Combine both plots
  p5 = ggarrange(pgs, pgsM, ncol = 2, nrow = 1, legend = "none")
  p5
  pdf(file = paste("plots_by_norm_grouped_mean-by-chrom/", aname,sep = ""), width =12, height = 12)
  
  plot(p5)
  dev.off()
  
}

#######Group by individuals and chromosome - normalised by cov

pop <- read.table("../bowtie_practise.txt") #list of sample names
slist <- pop$V1
str(slist)

for (sam in slist){
  name <- paste(sam,"_AchrnegJOIN.csv", sep = "")
  print (name) 
  
  AchrnegJOIN <- read.csv(name, header = TRUE)
  
  #Add columns with the names
  AchrnegJOIN$samples <- sam
  AchrnegJOIN[c('Number', 'Sample_name')] <- str_split_fixed(AchrnegJOIN$samples, "_", 2)
  
  if (sam == "S0001_G50516H"){
  AchrnegJOIN_all <- AchrnegJOIN
  }
  
  else {
    AchrnegJOIN_all <- rbind(AchrnegJOIN, AchrnegJOIN_all)
  }

}

write.csv(Mchr_all, file = "Mchr_all.csv")
write.csv(Mchr_all, file = "AchrnegJOIN_all.csv")


p <- ggplot(Cluster, aes(Sample_name, Mean_cov))
p1 <- p + geom_line(col = "grey60")
p2 <- p1 + geom_point(aes(colour=genome))
p3 <- p2 + theme(axis.text.x=element_text(angle = 90, hjust = 1, size=6)) + xlab("Chromosome") + ylab("mean percentage coverage mapped to each chromosome")
p4 <- p3 +  scale_color_manual(values=mycolors_change)
p5 <- p4 + facet_wrap(~Chromosome)
p5

p <- ggplot() + geom_bar(data=Mchr01,aes(x=start,y=covNormByMean), stat = "identity",color="firebrick3")
pg <- p + geom_bar(data=Achr01negJOIN,aes(x=inM_start,y=covbynorm_neg), stat = "identity",color="dodgerblue")
pgx <- pg + coord_cartesian(ylim = c(-3,3)) + theme(axis.text = element_text(size = 5))  
pgxl <- pgx + geom_hline(yintercept=0, linetype="solid", color = "gray40")
pgxk <- pgxl + geom_hline(yintercept=maMn, linetype="dashed", color = "gray30", alpha=0.3)
pgxm <- pgxk + geom_hline(yintercept=-maAn, linetype="dashed", color = "gray30", alpha=0.3)
pgsw <- pgxm + facet_wrap(~facet, ncol=1, strip.position = "right")
pgsM <- pgsw + theme_classic()
pgsM


###-----
##Plotting the difference between the normalised mean by chromosome
###Below doesn't work

for (sam in slist) {
  name <- paste(sam,"_unique.bam.bedgraph_mean_cov.txt",sep = "")
  print (name) 
  X <- sam
  #all <- read.delim ("S0001_G50516H_pp.bam_rep.bedgraph_median_cov.txt", header = FALSE)
  all<- read.delim (name, header = FALSE)
  colnames(all) <- c("chrom","start","end","cov")
  head(all)
  
  # to set window with no coverage of 0, these are input as "." 
  all$cov <- as.numeric(all$cov) 
  all$cov[is.na(all$cov)] <- 0
  
  print (sam)
  #colnames(all) <- c("chrom","start","end","cov")
  
  cov_chr <- aggregate(all[,4], list(all$chrom), mean)
  colnames(cov_chr) <- c("chrom", "mean_cov")
  
  all <- merge(all, cov_chr, by = "chrom")
  
  ccname <- paste(sam,"_chr_cov.csv",sep = "")
  write.csv(cov_chr,ccname)
  
  ma <- mean(all$cov)
  print(ma)
  mean_list <- append(mean_list,c(ma))
  
  all$covNormByMean <- all$cov/all$mean_cov ###relative coverage
  check <- mean(all$covNormByMean)
  print(check) #Should =1
  
  # get A genome chromosomes
  Achr01<- all %>% filter(grepl('^ChrA', chrom))
  maA <-  mean(Achr01$cov)
  normA_list <- append(normA_list,c(maA))
  maAn <- mean(Achr01$covNormByMean)
  Amean_list <- append(Amean_list,c(maAn))
  
  # get M genome chromosomes
  Mchr01<- all %>% filter(grepl('^ChrM', chrom))
  maM <- mean(Mchr01$cov)
  normM_list <- append(normM_list,c(maM))
  maMn <- mean(Mchr01$covNormByMean)
  Mmean_list <- append(Mmean_list,c(maMn))
  
  ##Define function 
  all.neg <- function(x) -1*abs(x)
  
  ## make B negative and remove B from chromosome -
  
  Mchr01neg <- cbind(Mchr01,"covbynorm_neg" = all.neg(Mchr01$covNormByMean))
  #Mchr01neg <- cbind(Mchr01,"cov_neg" = all.neg(Mchr01$cov))
  
  Achr01neg <- cbind(Achr01,"covbynorm_neg" = all.neg(Achr01$covNormByMean))
  #Achr01neg <- cbind(Achr01,"cov_neg" = all.neg(Achr01$cov))
  
  Mnes <- Mchr01neg %>%
    mutate(chrom = str_replace(chrom,"^ChrM","Chr"))
  
  Mpos <- Mchr01 %>%
    mutate(chrom = str_replace(chrom,"^ChrM","Chr"))
  
  Apos <- Achr01 %>% 
    mutate(chrom = str_replace(chrom, "^ChrA", "Chr"))
  
  Anes <- Achr01neg %>% 
    mutate(chrom = str_replace(chrom, "^ChrA", "Chr"))
  
  #print (head(Mnes))
  #print (head(Mchr01))
  
  #######################
  ###replot with syntheny
  #######################
  #read minimap alignments: 
  minimap <- read.delim("../../M_100kbwindows-over-A.longest_alignmentotal.paf", header=F)
  minimap <- minimap[,c(1,6,8,9)]
  colnames(minimap) <- c("uniqID","inAA_chr","inAA_start","inAA_stop")
  minimap <- minimap[- grep("scaffoldA", minimap$inAA_chr),]
  #create id field in Mchr01neg
  Mchr01neg$uniqID <- paste0(Mchr01neg$chrom,":",Mchr01neg$start,"-",Mchr01neg$end)
  
  #join left both dataframes
  Mchr01negJOIN <- merge(x=Mchr01neg, y=minimap, by="uniqID", all.x = FALSE) #false; so as to ignore Ms without A homologous
  
  Achr01$facet <- Achr01$chrom
  Mchr01negJOIN$facet <- Mchr01negJOIN$inAA_chr
  
  aname <- paste(sam,"cov_byA.pdf",sep = "")
  #aname <- 'file_cov_byA'
  
  print (aname)
  print (head(Achr01))

  check <- Achr01$covNormByMean - Mchr01negJOIN$covbynorm_neg
  
 #}
    
  ####Cov by mean - normalised
  p <- ggplot() + geom_bar(data=Achr01,aes(x=start,y=covNormByMean), stat = "identity",color="dodgerblue")
  pg <- p + geom_bar(data=Mchr01negJOIN,aes(x=inAA_start,y=covbynorm_neg), stat = "identity",color="firebrick3")
  pgx <- pg + coord_cartesian(ylim = c(-2.5,2.5)) + theme(axis.text = element_text(size = 5))  
  pgxl <- pgx + geom_hline(yintercept=0, linetype="solid", color = "gray40")
  pgxk <- pgxl + geom_hline(yintercept=maAn, linetype="dashed", color = "gray30", alpha=0.3)
  pgxm <- pgxk + geom_hline(yintercept=-maMn, linetype="dashed", color = "gray30", alpha=0.3)
  pgsw <- pgxm + facet_wrap(~facet, ncol=1, strip.position = "right")
  pgs <- pgsw + theme_classic()
  pgs
  
  pdf(file = paste("plots_by_norm_mean-by-chrom/", aname,sep = ""), width =8, height = 12)
  plot(pgs)
  dev.off()
  
  write.csv(Achr01, file = paste(X, "_Achr.csv", sep = ""))
  #write.csv(Mchr01negJOIN, file = paste(X, "_MchrnegJOIN.csv", sep = ""))
  
  ###SAME WITH M AS REFERENCE (ALSO SWAPPING TOP AND BOTTOM BUT NOT COLOURS):
  minimapM <- read.delim("../../A_100kbwindows-over-M.longest_alignmentotal.paf", header=F)
  minimapM <- minimapM[,c(1,6,8,9)]
  colnames(minimapM) <- c("uniqID","inM_chr","inM_start","inM_stop")
  minimapM <- minimapM[- grep("scaffoldM", minimapM$inM_chr),]
  #create id field in Bchr01neg
  Achr01neg$uniqID <- paste0(Achr01neg$chrom,":",Achr01neg$start,"-",Achr01neg$end)
  
  #join left both dataframes
  Achr01negJOIN <- merge(x=Achr01neg, y=minimapM, by="uniqID", all.x = FALSE) #false; so as to ignore As without M homologous
  
  Achr01negJOIN$facet <- Achr01negJOIN$inM_chr
  Mchr01$facet <- Mchr01$chrom
  
  bname <- paste(sam,"cov_byM.pdf",sep = "")
  #bname <- 'file_cov_byM'
  
  ####Cov by mean - normalised
  p <- ggplot() + geom_bar(data=Mchr01,aes(x=start,y=covNormByMean), stat = "identity",color="firebrick3")
  pg <- p + geom_bar(data=Achr01negJOIN,aes(x=inM_start,y=covbynorm_neg), stat = "identity",color="dodgerblue")
  pgx <- pg + coord_cartesian(ylim = c(-2.5,2.5)) + theme(axis.text = element_text(size = 5))  
  pgxl <- pgx + geom_hline(yintercept=0, linetype="solid", color = "gray40")
  pgxk <- pgxl + geom_hline(yintercept=maMn, linetype="dashed", color = "gray30", alpha=0.3)
  pgxm <- pgxk + geom_hline(yintercept=-maAn, linetype="dashed", color = "gray30", alpha=0.3)
  pgsw <- pgxm + facet_wrap(~facet, ncol=1, strip.position = "right")
  pgsM <- pgsw + theme_classic()
  pgsM
  
  pdf(file = paste("plots_by_norm_mean-by-chrom/", bname,sep = ""), width =8, height = 12)
  plot(pgsM)
  dev.off()
  
  write.csv(Mchr01, file = paste(X, "_Mchr.csv", sep = ""))
  #write.csv(Achr01negJOIN, file = paste(X, "_AchrnegJOIN.csv", sep = ""))
  
  ###Combine both plots
  p5 = ggarrange(pgs, pgsM, ncol = 2, nrow = 1, legend = "none")
  p5
  pdf(file = paste("plots_by_norm_grouped_mean-by-chrom/", aname,sep = ""), width =12, height = 12)
  
  plot(p5)
  dev.off()
  
}
