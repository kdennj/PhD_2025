library(readxl) 
library(qqman)
library(ggplot2)
library(dplyr)
library(grDevices)
library(knitr)
library(ggpubr)
library(RColorBrewer)

### Whole determinate files
det_gwas_Res_F <- read.csv("../../pca3_thin5_bytrait/Determinate/GAPIT.Association.GWAS_Results.FarmCPU.Determinate.csv")
det_gwas_Res_M <- read.csv("../../pca3_thin5_bytrait/Determinate/GAPIT.Association.GWAS_Results.MLM.Determinate.csv")
det_gwas_Res_B <- read.csv("../../pca3_thin5_bytrait/Determinate/GAPIT.Association.GWAS_Results.Blink.Determinate.csv")

#Add group column
det_gwas_Res_B$Model='Blink'
det_gwas_Res_F$Model='FarmCPU'
det_gwas_Res_M$Model='MLM'

#Combine dataframes
gwasResults <- rbind(det_gwas_Res_B, det_gwas_Res_F, det_gwas_Res_M)

##### Whole panel Manhattan plot -----------------------------------
x <- det_gwas_Res_B 
x <- det_gwas_Res_F
x <- det_gwas_Res_M

#Repeat for all x
don_B_W <- x %>% 
  
  # Compute chromosome size
  group_by(Chr) %>% 
  summarise(chr_len=max(Pos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(x, ., by=c("Chr"="Chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chr, Pos) %>%
  mutate( BPcum=Pos+tot)

axisdf_B_W = don_B_W %>%
  group_by(Chr) %>%
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

gwasResults_W <- rbind(don_B_W, don_F_W, don_M_W)

gwasResults_W$Sizes <- -log10(gwasResults_W$P.value)
gwasResults_W$size <- scales::rescale(gwasResults_W$Sizes, to = c(1, 7))

sig_MFB_W <- -log10(0.05 / nrow(gwasResults_W))

##Add vertical lines for matching with Andean or Whole panel
new <- gwasResults_W[gwasResults_W$Chr == "Chr01", ] 
new2 <- new[new$Pos > 6531240,]

#Determinate A matching with Andean or Whole panel
Det_W_MFB <- c(459848492, 398660040, 362380452, 311617117, 44941717, 42413434, 6531249)

## Plotting
pW <- ggplot(gwasResults_W, aes(x=BPcum, y=-log10(P.value), shape = Model)) +
  
  # Show all points
  geom_point( aes(color=as.factor(Chr)), alpha=0.8, size=gwasResults_W$size, stroke = 1.5) +
  #scale_color_manual(values = rep(c("#ff6666", "#ffbd55", "#9de24f", "#87cefa", "#cf9dff"), 11 )) +
  scale_color_manual(values = rep(c("#87cefa", "darkgrey"), 11 )) +
  
  #Bonferonni cutoff
  geom_hline(
    yintercept = sig_MFB_W, color = "forestgreen", linewidth = 1.4) +
  
  ##Add vertical lines for matching with Andean or Whole panel  
  geom_vline(xintercept = Det_W_MFB, color = "black", linetype = "dotted", linewidth = 1.4) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$Chr, breaks= axisdf$center, expand = c(0, 0) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,55) ) +     # remove space between plot area and x axis
  
  guides(colour = "none") +
  
  # Custom the theme:
  theme_classic() +
  theme( 
    axis.title.x=element_blank(),
    #legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size=18),
    axis.title.y = element_text(size=18),
    legend.text = element_text(size=16),
    legend.title = element_text(size=18)
  ) +
  guides(shape = guide_legend(override.aes = list(size = 5)))
pW

png("det_gwas_MFB_manhattan_new.png",  width = 1500, height = 600)
plot(pW)
dev.off()