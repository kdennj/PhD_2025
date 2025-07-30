library("RColorBrewer")
library("ggplot2")

library(dplyr)
library(tidyr)
library(stringr)

rm(list = ls())

#Following what Janet used - median
setwd("~/OneDrive/STATS/22.03 Chromosome coverage/bwa/23.01-median_cov")

setwd("STATS/22.03 Chromosome coverage/bwa/23.01-mean_cov")

# script to read in coverage files  - median coverage over 100,000bp windows for A adn B genome
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
  name <- paste(sam,"_pp.bam_rep.bedgraph_mean_cov_rep.txt",sep = "")
  print (name) 
  #all <- read.delim ("S0001_G50516H_pp.bam_rep.bedgraph_median_cov.txt", header = FALSE)
  all<- read.delim (name, header = FALSE)
  colnames(all) <- c("chrom","start","end","cov")
  head(all)

  # to set window with no coverage of 0, these are input as "." 
  all$cov <- as.numeric(all$cov) 
  all$cov[is.na(all$cov)] <- 0
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
  
  pdf(file = paste("Plots_by_cov/", bname,sep = ""),width =8, height = 10)
  plot(pgsM)
  dev.off()

  ###Combine both plots
  p5 = ggarrange(pgs, pgsM, ncol = 2, nrow = 1, legend = "none")
  p5
  pdf(file = paste("plots_by_cov_grouped/", aname,sep = ""), width =12, height = 12)
  
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

#Following what Janet used - median
setwd("/STATS/22.03 Chromosome coverage/bwa/23.01-median_cov")

setwd("STATS/22.03 Chromosome coverage/bwa/23.01-mean_cov")

rm(list = ls())

# script to read in coverage files  - median coverage over 100,000bp windows for A adn B genome
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
  name <- paste(sam,"_pp.bam_rep.bedgraph_mean_cov_rep.txt",sep = "")
  print (name) 
  #all <- read.delim ("S0001_G50516H_pp.bam_rep.bedgraph_median_cov.txt", header = FALSE)
  all<- read.delim (name, header = FALSE)
  colnames(all) <- c("chrom","start","end","cov")
  head(all)
  
  # to set window with no coverage of 0, these are input as "." 
  all$cov <- as.numeric(all$cov) 
  all$cov[is.na(all$cov)] <- 0
  
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
  
  pdf(file = paste("plots_by_norm/", aname,sep = ""), width =8, height = 12)
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
  
  pdf(file = paste("plots_by_norm/", bname,sep = ""), width =8, height = 12)
  plot(pgsM)
  dev.off()
  
  
  ###Combine both plots
  p5 = ggarrange(pgs, pgsM, ncol = 2, nrow = 1, legend = "none")
  p5
  pdf(file = paste("plots_by_norm_grouped_together/", aname,sep = ""), width =12, height = 12)
  
  plot(p5)
  dev.off()
  
}



