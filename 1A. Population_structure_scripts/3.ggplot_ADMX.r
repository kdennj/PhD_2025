library(grid)
library(ggplot2)
library(tidyr)
library("RColorBrewer")
library(plyr)
library(readxl)
library(dplyr)
library(forcats)
library("gridExtra")
library("cowplot")
library(patchwork)
library(ggpubr)


pcaPalette <- c("#2121D9","#9999FF","#FFFB23","#FF9326","#A945FF","#0089B2")
safe_pal_8 <- c("#000000", "#888888", "#5EC962", "#117733", "#332288", "#6699CC", "#CC5500","#AA4499")

rm(list = ls())

######### plot K6 ########
K6_clustering <- read_excel("K6_Q_files/K6_clustering_PCA-phenotypes over 0.7 ADMX.xlsx")
K6_clustering$Sample_name <- factor(K6_clustering$Sample_name, levels = K6_clustering$Sample_name)

#### Re-order the files for plotting 
K6_clustering_order <- K6_clustering[order(K6_clustering$New_order),]

##Change names of columns
colnames(K6_clustering_order) <- c('Sample_name', 'A1', 'C1', 'M1', 'C-EP', 'M2', 'C2', 'Admixed_K6_C', 'New_order')
K6rg <- gather(K6_clustering_order, key, value, -Sample_name, -New_order, -Admixed_K6_C)

K6rg$Sample_name <- as.character(K6rg$Sample_name)
K6rg$Sample_name <- factor(K6rg$Sample_name, levels = unique(K6rg$Sample_name))

match_PCA <- c("#117733" , "#332288" , "#5EC962" , "#AA4499" , "#6699CC" ,"#CC5500") 

####### plot pop proportion ############
cog <- ggplot(data=K6rg, aes(x=Sample_name, y=value, fill=key)) + 
  geom_bar(position="stack", stat="identity", show.legend = F)

cog3 <- cog +theme(axis.text.x=element_text(angle = 55, hjust = 1, size=4)) +
  xlab ("Genotype") + ylab("Q value") + scale_y_continuous(expand = c(0,0)) 

cog4 <- cog3 + scale_fill_manual(values=match_PCA) #+ geom_vline(xintercept = int)
        #  xintercept = K6rg$Admixed_K6_C)

cog5 <- cog4 + facet_grid(~fct_relevel(Admixed_K6_C, 'C-EP', 'A1', 'Admx_A', 'C1', 'C2', 'Admx_AM', 'M1', 'Admx_M', 'M2'), scales = "free_x", space = "free")
cog6 <- cog5 + theme(panel.spacing = unit(0.2, "lines"), strip.text.x = element_text(size = 6))

pdf(file = "ggplot_figures/plot_clusterK6_ggplot_legend_facetgrid_wout_legend3.pdf",width=12,height=2.5)
cog6
dev.off()

##########Plot for K2# 
K2_clustering <- read_excel("K2_Q_files/K2_clustering_PCA-phenotypes over 0.7 ADMX.xlsx")

K2_clustering$Sample_name <- factor(K2_clustering$Sample_name, levels = K2_clustering$Sample_name)

#### Re-order the files for plotting 
K2_clustering_order <- K2_clustering[order(K2_clustering$New_order),]

#Change column names
colnames(K2_clustering_order) <- c('Sample_name', 'Admixed_K6_C', 'Mesoamerican', 'Andean', 'New_order')

#K2rg <- gather(K2_clustering_order, key, value, -Sample_name, -Order_ADMX_K2, -Admixed_K6_C)
K2rg <- gather(K2_clustering_order, key, value, -Sample_name, -New_order, -Admixed_K6_C)

match_PCA <- c("#000000" , "#D3D3D3") 

K2rg$Sample_name <- as.character(K2rg$Sample_name)
K2rg$Sample_name <- factor(K2rg$Sample_name, levels = unique(K2rg$Sample_name))

p <- ggplot(data=K2rg, aes(x=Sample_name, y=value, fill=key)) + 
  geom_bar(stat="identity", show.legend = F)

p2 <- p +theme(axis.text.x=element_text(angle = 55, hjust = 1, size=4)) +
  xlab ("Genotype") + ylab("Q value") + scale_y_continuous(expand = c(0,0)) 

p3 <- p2 + scale_fill_manual(values=match_PCA) #+
  #geom_vline(xintercept = int, colour = "white") #  xintercept = K6rg$Admixed_K6_C)

p4 <- p3 + facet_grid(~fct_relevel(Admixed_K6_C, 'C-EP', 'A1', 'Admx_A', 'C1', 'C2', 'Admx_AM', 'M1', 'Admx_M', 'M2'), scales = "free_x", space = "free")

#Edited for gggarrange
p5 <- p4 + theme(panel.spacing = unit(0.2, "lines"), strip.text.x = element_text(size = 6), axis.title.x=element_blank(),
                 axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.background = element_rect(colour="white", fill="white")) 


pdf(file = "ggplot_figures/plot_clusterK2_ggplot_facetwrap_blank_xaxis.pdf",width=10,height=2)
p5
dev.off()

