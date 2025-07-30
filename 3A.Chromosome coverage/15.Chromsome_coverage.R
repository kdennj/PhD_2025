library("RColorBrewer")
library("ggplot2")

library(dplyr)
library(tidyr)
library(stringr)

rm(list = ls())

#### plot median
setwd("STATS/22.03 Chromosome coverage/median_cov")
getwd()


all <- read.delim ("S0001_G50516H_pp.bam.bedgraph_median_cov.txt", header = FALSE)
colnames(all) <- c("chrom","start","end","cov")

#Janet divided cov/per/sam by number of samples
#allc <- cbind(practise,"covPerSam" = practise$cov/22)

#First time didn't normalise for coverage 
practise <- read.delim ("S0001_G50516H_pp.bam.bedgraph_median_cov.txt", header = FALSE)
colnames(practise) <- c("chrom","start","end","cov") 

practise$cov[is.na(practise$cov)] <- 0
practise$cov <- as.numeric(practise$cov)

Achr01<- practise %>% filter(grepl('ChrA', chrom))
mean(Achr01$cov)
#[1] 9.647299

# get M genome chromosomes
Mchr01 <- practise %>% filter(grepl('ChrM', chrom))
mean(Mchr01$cov)
#[1] 8.502892

## make B negative and remove B from chromosome
all.neg <- function(x) -1*abs(x)

Mchr01neg <- cbind(Mchr01,"cov_neg" = all.neg(Mchr01$cov))

Mnes <- Mchr01neg %>%
  mutate(chrom = str_replace(chrom,"ChrM","Chr"))

Apos <- Achr01 %>% 
  mutate(chrom = str_replace(chrom, "ChrA", "Chr"))


## now do both A and B
p <- ggplot() + geom_bar(data=Apos,aes(x=start,y=cov), stat = "identity",color="dodgerblue")
pg <- p + geom_bar(data=Mnes,aes(x=start,y=cov_neg), stat = "identity",color="firebrick3")
pgx <- pg + coord_cartesian(ylim = c(-15, 15)) + theme(axis.text = element_text(size = 5))    
pgxl <- pgx + geom_hline(yintercept=0, linetype="solid", color = "gray40")
pgxk <- pgxl + geom_hline(yintercept=9.65, linetype="dashed", color = "gray30")
pgxm <- pgxk + geom_hline(yintercept=-8.50, linetype="dashed", color = "gray30")
pgsw <- pgxm + facet_wrap(~chrom,ncol=1,strip.position = "right")
pgs <- pgsw + theme_classic()

pdf(file = "1_AMchr_100kbcov.pdf",width =8, height = 12)
pgs
dev.off()


########################## Normalise coverage #######################################

# initialise empty lists
mean_list <- c()
normA_list <- c()
normM_list <- c()
Amean_list <- c()
Mmean_list <- c()

all$cov[is.na(all$cov)] <- 0
all$cov <- as.numeric(all$cov)
cov_chr <- all %>% group_by(chrom) %>%
  summarize(cov = mean(cov))
write.csv(cov_chr)
ma <- mean(all$cov)
print(ma)
#[1] 9.288651

all$covNormByMean <- all$cov/mean(all$cov)  ###relative coverage
check <- mean(all$covNormByMean)
mean(all$covNormByMean)
print(check)  #Should =1

# get A genome chromosomes
Achr01<- all %>% filter(grepl('ChrA', chrom))
maA <-  mean(Achr01$cov)
maAn <- mean(Achr01$covNormByMean)

print(maA)
#[1] 9.647299

print(maAn)
#[1] 1.038611

# get M genome chromosomes
Mchr01<- all %>% filter(grepl('ChrM', chrom))
maM <- mean(Mchr01$cov)
maMn <- mean(Mchr01$covNormByMean)

print(maM)
#[1] 8.502892

print(maMn)
#[1] 0.9154065

## make M negative and remove M and A from chromosome
all.neg <- function(x) -1*abs(x)

Mchr01neg <- cbind(Mchr01,"covbymean_neg" = all.neg(Mchr01$covNormByMean))

Mnes <- Mchr01neg %>%
  mutate(chrom = str_replace(chrom,"ChrM","Chr"))

Apos <- Achr01 %>% 
  mutate(chrom = str_replace(chrom, "ChrA", "Chr"))


print (head(Mnes))
print (head(Mchr01))


## now do both A and B
p <- ggplot() + geom_bar(data=Apos,aes(x=start,y=covNormByMean), stat = "identity",color="dodgerblue")
pg <- p + geom_bar(data=Mnes,aes(x=start,y=cov_neg), stat = "identity",color="firebrick3")
pgx <- pg + coord_cartesian(ylim = c(-5, 5)) + theme(axis.text = element_text(size = 5))    
pgxl <- pgx + geom_hline(yintercept=0, linetype="solid", color = "gray40")
pgxk <- pgxl + geom_hline(yintercept=1.04, linetype="dashed", color = "gray30")
pgxm <- pgxk + geom_hline(yintercept=-0.92, linetype="dashed", color = "gray30")
pgsw <- pgxm + facet_wrap(~chrom,ncol=1,strip.position = "right")
pgs <- pgsw + theme_classic()

pdf(file = "1_AMchr_100kbcov_normalised_5-5.pdf",width =8, height = 12)
pgs
dev.off()







all <- read.delim ("1_G50516H_pp.bam.bedgraph_median_cov.txt", header = FALSE)

colnames(all) <- c("chrom","start","end","cov")

aver_practise_cov <- all %>%                                        
  group_by(chrom) %>%                         
  summarise_at(vars(cov),              
               list(name = mean))

write.csv(aver_practise_cov,"covPer_chrom.csv")

aver_pp_cov <- all %>%                                        
  group_by(chrom) %>%                         
  summarise_at(vars(covPerSam),              
               list(name = mean))

write.csv(aver_pp_cov,"covPer_chrom.csv")


#######################
###replot with syntheny
#######################
#read minimap alignments:
minimap <- read.delim("M_100kbwindows-over-A.longest_alignmentotal.paf", header=F)
minimap <- minimap[,c(1,6,8,9)]
colnames(minimap) <- c("uniqID","inAA_chr","inAA_start","inAA_stop")
minimap <- minimap[- grep("scaffoldA", minimap$inAA_chr),]
#create id field in Mchr01neg
Mchr01neg$uniqID <- paste0(Mchr01neg$chrom,":",Mchr01neg$start,"-",Mchr01neg$end)

#join left both dataframes
Mchr01negJOIN <- merge(x=Mchr01neg, y=minimap, by="uniqID", all.x = FALSE) #false to ignore BBs without AA homologous

Achr01$facet <- Achr01$chrom
Mchr01negJOIN$facet <- Mchr01negJOIN$inAA_chr

aname <- paste(sam,"cov_byA.pdf",sep = "")

print (aname)
print (head(Achr01))


p <- ggplot() + geom_bar(data=Achr01,aes(x=start,y=covNormByMean), stat = "identity",color="dodgerblue")
pg <- p + geom_bar(data=Mchr01negJOIN,aes(x=inAA_start,y=cov_neg), stat = "identity",color="firebrick3")
pgx <- pg + coord_cartesian(ylim = c(-2,2)) + theme(axis.text = element_text(size = 5))  
pgxl <- pgx + geom_hline(yintercept=0, linetype="solid", color = "gray40")
pgxk <- pgxl + geom_hline(yintercept=1.04, linetype="dashed", color = "gray30", alpha=0.3)
pgxm <- pgxk + geom_hline(yintercept=-0.92, linetype="dashed", color = "gray30", alpha=0.3)
pgsw <- pgxm + facet_wrap(~facet, ncol=1, strip.position = "right")
pgs <- pgsw + theme_classic()
pgs

pdf(file = aname,width =8, height = 12)
plot(pgs)
dev.off()


###SAME WITH M AS REFERENCE (ALSO SWAPPING TOP AND BOTTOM BUT NOT COLOURS):

Achr01neg <- cbind(Achr01,"covbymean_neg" = all.neg(Achr01$covNormByMean))

minimapM <- read.delim("A_100kbwindows-over-M.longest_alignmentotal.paf", header=F)
minimapM <- minimapM[,c(1,6,8,9)]
colnames(minimapM) <- c("uniqID","inM_chr","inM_start","inM_stop")
minimapM <- minimapM[- grep("scaffoldM", minimapM$inM_chr),]

minimapMChr03 <- minimapM[grep("ChrM03", minimapM$inM_chr),]
minimapMChr02 <- minimapM[grep("ChrM02", minimapM$inM_chr),]


#create id field in Bchr01neg
Achr01neg$uniqID <- paste0(Achr01neg$chrom,":",Achr01neg$start,"-",Achr01neg$end)

#join left both dataframes
Achr01negJOIN <- merge(x=Achr01neg, y=minimapM, by="uniqID", all.x = FALSE) #false to ignore BBs without AA homologous
Achr01negJOIN$facet <- Achr01negJOIN$inM_chr

Mchr01$facet <- Mchr01$chrom

bname <- paste(sam,"cov_byB.pdf",sep = "")

p <- ggplot() + geom_bar(data=Mchr01,aes(x=start,y=covNormByMean), stat = "identity",color="firebrick3")
pg <- p + geom_bar(data=Achr01negJOIN,aes(x=inM_start,y=cov_neg), stat = "identity",color="dodgerblue")
pgx <- pg + coord_cartesian(ylim = c(-2,2)) + theme(axis.text = element_text(size = 5))  
pgxl <- pgx + geom_hline(yintercept=0, linetype="solid", color = "gray40")
pgxk <- pgxl + geom_hline(yintercept=1.04, linetype="dashed", color = "gray30", alpha=0.3)
pgxm <- pgxk + geom_hline(yintercept=-0.92, linetype="dashed", color = "gray30", alpha=0.3)
pgsw <- pgxm + facet_wrap(~facet,ncol=1,strip.position = "right")
pgs <- pgsw + theme_classic()
pgs

pdf(file = bname,width =8, height = 12)
plot(pgs)
dev.off()


#######################
###replot with syntheny with coverage (not coverage by mean)
#######################
#read minimap alignments:
minimap <- read.delim("M_100kbwindows-over-A.longest_alignmentotal.paf", header=F)
minimap <- minimap[,c(1,6,8,9)]
colnames(minimap) <- c("uniqID","inAA_chr","inAA_start","inAA_stop")
minimap <- minimap[- grep("scaffoldA", minimap$inAA_chr),]
#create id field in Mchr01neg
Mchr01neg$uniqID <- paste0(Mchr01neg$chrom,":",Mchr01neg$start,"-",Mchr01neg$end)

#join left both dataframes
Mchr01negJOIN <- merge(x=Mchr01neg, y=minimap, by="uniqID", all.x = FALSE) #false to ignore BBs without AA homologous

Achr01$facet <- Achr01$chrom
Mchr01negJOIN$facet <- Mchr01negJOIN$inAA_chr

aname <- paste(sam,"cov_byA.pdf",sep = "")

print (aname)
print (head(Achr01))


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


###SAME WITH M AS REFERENCE (ALSO SWAPPING TOP AND BOTTOM BUT NOT COLOURS):

Achr01neg <- cbind(Achr01,"cov_neg" = all.neg(Achr01$covNormByMean))

minimapM <- read.delim("A_100kbwindows-over-M.longest_alignmentotal.paf", header=F)
minimapM <- minimapM[,c(1,6,8,9)]
colnames(minimapM) <- c("uniqID","inM_chr","inM_start","inM_stop")
minimapM <- minimapM[- grep("scaffoldM", minimapM$inM_chr),]

minimapMChr03 <- minimapM[grep("ChrM03", minimapM$inM_chr),]
minimapMChr02 <- minimapM[grep("ChrM02", minimapM$inM_chr),]


#create id field in Bchr01neg
Achr01neg$uniqID <- paste0(Achr01neg$chrom,":",Achr01neg$start,"-",Achr01neg$end)

#join left both dataframes
Achr01negJOIN <- merge(x=Achr01neg, y=minimapM, by="uniqID", all.x = FALSE) #false to ignore BBs without AA homologous
Achr01negJOIN$facet <- Achr01negJOIN$inM_chr

Mchr01$facet <- Mchr01$chrom

bname <- paste(sam,"cov_byB.pdf",sep = "")

p <- ggplot() + geom_bar(data=Mchr01,aes(x=start,y=covNormByMean), stat = "identity",color="firebrick3")
pg <- p + geom_bar(data=Achr01negJOIN,aes(x=inM_start,y=cov_neg), stat = "identity",color="dodgerblue")
pgx <- pg + coord_cartesian(ylim = c(-2,2)) + theme(axis.text = element_text(size = 5))  
pgxl <- pgx + geom_hline(yintercept=0, linetype="solid", color = "gray40")
pgxk <- pgxl + geom_hline(yintercept=1.04, linetype="dashed", color = "gray30", alpha=0.3)
pgxm <- pgxk + geom_hline(yintercept=-0.92, linetype="dashed", color = "gray30", alpha=0.3)
pgsw <- pgxm + facet_wrap(~facet,ncol=1,strip.position = "right")
pgs <- pgsw + theme_classic()
pgs

pdf(file = bname,width =8, height = 12)
plot(pgs)
dev.off()






















