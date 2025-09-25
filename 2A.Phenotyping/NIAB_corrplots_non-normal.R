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
library("corrplot")
library(ggplot2)
library(gridExtra)
library(grid)

setwd("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Bean project/23.03-Field-trial-NIAB/Licor/")

### Corrplot non-normal distribution -----------------
###This is for non-normal corrplots 
harvest <- read_excel("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Bean project/23.03-Field-trial-NIAB/KDJ Bean panel randomisation.xlsx", sheet = "harvest_edit")
all_df <- read.xlsx("Plots/corrplots/all_dates_licor.xlsx")
seeds <- read.xlsx("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Bean project/Common bean samples/22.02 Seed phenotypes/22.02-updated-phenotypes-R.xlsx")
seeds <- seeds[, c(1,6, 20)]
seeds$Sample_name <- paste(seeds$Sample_name, "D", sep = "_")
colnames(seeds)[1] <- c("Accession_T")
colnames(seeds)[2] <- c("Major_seed_colour")

harvest$Accession_T <- paste(harvest$Accession, harvest$Treat, sep = "_") 
harvest <- harvest[, c(22, 6, 12:14)]
#harvest$GH <- gsub('Determinate_bush', '1', harvest$GH)
#harvest$GH <- gsub('Indeterminate_bush', '2', harvest$GH)
#harvest$GH <- gsub('Indeterminate_climbing', '3', harvest$GH)
harvest[, -c(1, 2)] <- lapply(harvest[, -c(1, 2)], as.numeric)
harvest <- aggregate(. ~ Accession_T + GH, 
                     data = harvest, 
                     FUN = mean, 
                     na.rm = TRUE)

#Need for popstruc information
date28.07 <- read.xlsx("28.07.23.merged.xlsx")
#which(colnames(date28.07) == "Admixed_K6_C")
date28.07$Accession_T <- paste(date28.07$Accession, date28.07$Treat, sep = "_")
popstruc <- date28.07[, c(111, 108)]
# popstruc$Admixed_K6_C <- gsub('A1', '1', popstruc$Admixed_K6_C)
# popstruc$Admixed_K6_C <- gsub('Admx_A', '2', popstruc$Admixed_K6_C)
# popstruc$Admixed_K6_C <- gsub('Admx_AM', '3', popstruc$Admixed_K6_C)
# popstruc$Admixed_K6_C <- gsub('Admx_M', '4', popstruc$Admixed_K6_C)
# popstruc$Admixed_K6_C <- gsub('C-EP', '5', popstruc$Admixed_K6_C)
# popstruc$Admixed_K6_C <- gsub('C1', '6', popstruc$Admixed_K6_C)
# popstruc$Admixed_K6_C <- gsub('C2', '7', popstruc$Admixed_K6_C)
# popstruc$Admixed_K6_C <- gsub('M1', '8', popstruc$Admixed_K6_C)
# popstruc$Admixed_K6_C <- gsub('M2', '9', popstruc$Admixed_K6_C)
# popstruc$Admixed_K6_C <- as.numeric(popstruc$Admixed_K6_C)
# popstruc$Admixed_K6_C <- replace_na(popstruc$Admixed_K6_C, 0) 
# popstruc <- unique(popstruc)

strategies <- read.xlsx("mean_weights_pods_foliar_all_correct.xlsx")
strategies <- strategies[, c(1, 17, 19)]
colnames(strategies)[2] <- c("Strategy")
colnames(strategies)[3] <- c("Recovery")
# strategies$Strategies <- gsub('DB control growth', '1', strategies$Strategies)
# strategies$Strategies <- gsub('Drought susceptible', '2', strategies$Strategies)
# strategies$Strategies <- gsub('No classification', '3', strategies$Strategies)
# strategies$Strategies <- gsub('Prioritised yield', '4', strategies$Strategies)
# strategies$Strategies <- gsub('Saver', '5', strategies$Strategies)
# strategies$Strategies <- gsub('Spender', '6', strategies$Strategies)
# strategies$Strategies <- gsub('Stay-green', '7', strategies$Strategies)
# strategies$Strategies <- as.numeric(strategies$Strategies)
# strategies$Strategies <- replace_na(strategies$Strategies, 0) 

all_join <- join(all_df, popstruc, by = "Accession_T") %>%
  unique()

all_join <- join(all_join, harvest, by = "Accession_T") %>%
  unique()

all_join <- join(all_join, strategies, by = "Accession_T") %>%
  unique()

stay_green <- all_join[grepl("Saver", all_join$Column2),]
stay_green <- stay_green[, -c(17:25)]
stay_green <- aggregate(. ~ Accession_T + date + Column2, 
                    data = stay_green, 
                    FUN = mean, 
                    na.rm = TRUE)
stay_green <- stay_green[grepl("date10.08", stay_green$date),]
x <- stay_green$Admixed_K6_C %>% unique()
x

recovery <- all_join[grepl("Recovery", all_join$Descriptor_recovery),]
recovery <- recovery[, -c(17:26)]
recovery[, -c(1, 11, 17)] <- lapply(recovery[, -c(1, 11, 17)], as.numeric)
recovery <- aggregate(. ~ Accession_T + date + Descriptor_recovery, 
                        data = recovery, 
                        FUN = mean, 
                        na.rm = TRUE)
recovery <- recovery[grepl("date10.08", recovery$date),]
x <- recovery$Admixed_K6_C %>% count()
x

#all_join <- all_join[, c(-11)]

all_join[, -c(1, 11, 12, 13, 17, 18)] <- lapply(all_join[, -c(1, 11, 12, 13, 17, 18)], as.numeric) #11:13, 17

result <- aggregate(. ~ Accession_T + date + GH + Admixed_K6_C + Strategy + Recovery, #GH + Admixed_K6_C + Strategy
                    data = all_join, 
                    FUN = mean, 
                    na.rm = TRUE)

determinate <- result[grepl("1", result$GH),]

result <- result[!grepl("_C", result$Accession_T), ]
#Remove those which we haven't sequenced
accession_list <- c("G50516A_D", "G50516O_D", "G50516D_D", "G51285_D", "G51285H_D", "G8152_D", "G916_D", "G50997_D")  
result <- result[!result$Accession_T %in% accession_list, ]

result <- result[, c(-1, -2)]

?corrplot
M = cor(result, method = "spearman")
testRes = cor.mtest(M, conf.level = 0.95)

pdf("Plots/corrplots/non_normal/Drought_allsep_sequenced_corrplot_lower_hclust_weight_strat.pdf", height=7, width = 7)
corrplot(M, p.mat = testRes$p, method = 'color', diag = T, type = 'lower',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20', order = 'hclust', outline = T, 
         col = COL2('RdYlBu'), tl.col = 'black') #Col1 is for sequential colours, Col2 is for diverging
dev.off()

#With BBCH
date30.08 <- read.xlsx("30.08.23.merged.xlsx")
date23.08 <- read.xlsx("23.08.23.merged.xlsx")
date18.08 <- read.xlsx("18.08.23.merged.xlsx")
date10.08 <- read.xlsx("10.08.23.merged.xlsx")
date28.07 <- read.xlsx("28.07.23.merged.xlsx")

date30.08$Accession_T <- paste(date30.08$Accession, date30.08$Treat, sep = "_")
date23.08$Accession_T <- paste(date23.08$Accession, date23.08$Treat, sep = "_")
date18.08$Accession_T <- paste(date18.08$Accession, date18.08$Treat, sep = "_")
date10.08$Accession_T <- paste(date10.08$Accession, date10.08$Treat, sep = "_")
date28.07$Accession_T <- paste(date28.07$Accession, date28.07$Treat, sep = "_")

date30.08 <- date30.08[, c(111, 106)]
date23.08 <- date23.08[, c(111, 106)]
date18.08 <- date18.08[, c(111, 106)]
date10.08 <- date10.08[, c(111, 106)]
date28.07 <- date28.07[, c(111, 106)]
#which(colnames(date28.07) == "28/07/2023")

date30.08$date <- c("date30.08")
date23.08$date <- c("date23.08")
date18.08$date <- c("date18.08")
date10.08$date <- c("date10.08")
date28.07$date <- c("date28.07")

colnames(date30.08)[2] <- c("BBCH")
colnames(date23.08)[2] <- c("BBCH")
colnames(date18.08)[2] <- c("BBCH")
colnames(date10.08)[2] <- c("BBCH")
colnames(date28.07)[2] <- c("BBCH")

df_list <- list(date30.08, date23.08, date18.08, date10.08, date28.07)

BBCH <- bind_rows(df_list)

BBCH_res <- aggregate(. ~ Accession_T + date, 
                      data = BBCH, 
                      FUN = mean, 
                      na.rm = TRUE)

BBCH_res <- BBCH_res[!grepl("_C", BBCH_res$Accession_T), ]
BBCH_res <- BBCH_res[!BBCH_res$Accession_T %in% accession_list, ]

bbch_all_result <- join(BBCH_res, result, by = c("Accession_T", "date"))
#bbch_all_result <- bbch_all_result[, c(-1, -2)]

bbch_all_result <- join(bbch_all_result, seeds)
bbch_all_result <- join(bbch_all_result, strategies)

#Add country and GH
K6_clustering <- read_excel("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Bean project/STATS/22.01 Variant calling - GATK/22.11-ADMIXTURE-LD0.5/Q_files/K6_Q_files/K6_clustering_PCA-phenotypes over 0.7 ADMX.xlsx")
#which(colnames(K6_clustering) == 'photo sens')

more <- K6_clustering[, c(2, 9, 10, 47)]
colnames(more)[1] <- c("Accession")

bbch_all_result$Accession_T <- gsub("_D$", "", bbch_all_result$Accession_T)
colnames(bbch_all_result)[1] <- c("Accession")
bbch_all_result <- join(bbch_all_result, more, by = c("Accession"))

# bbch_all_result$Type <- gsub('Heirloom', '1', bbch_all_result$Type)
# bbch_all_result$Type <- gsub('Landrace', '2', bbch_all_result$Type)
# bbch_all_result$Type <- gsub('Wild', '3', bbch_all_result$Type)
# bbch_all_result$Country <- gsub('Argentina', '1', bbch_all_result$Country)
# bbch_all_result$Country <- gsub('Brazil', '2', bbch_all_result$Country)
# bbch_all_result$Country <- gsub('Chile', '3', bbch_all_result$Country)
# bbch_all_result$Country <- gsub('Colombia', '4', bbch_all_result$Country)
# bbch_all_result$Country <- gsub('Ecuador', '5', bbch_all_result$Country)
# bbch_all_result$Country <- gsub('Guatemala', '6', bbch_all_result$Country)
# bbch_all_result$Country <- gsub("Heirlooms\\(NA\\)", "7", bbch_all_result$Country)
# bbch_all_result$Country <- gsub('Honduras', '8', bbch_all_result$Country)
# bbch_all_result$Country <- gsub('Mexico', '9', bbch_all_result$Country)
# bbch_all_result$Country <- gsub('Peru', '10', bbch_all_result$Country)

# bbch_all_result$date <- gsub('date10.08', '1', bbch_all_result$date)
# bbch_all_result$date <- gsub('date18.08', '2', bbch_all_result$date)
# bbch_all_result$date <- gsub('date23.08', '3', bbch_all_result$date)
# bbch_all_result$date <- gsub('date30.08', '4', bbch_all_result$date)
# bbch_all_result$date <- gsub('date28.07', '5', bbch_all_result$date)


#bbch_all_result <- bbch_all_result[, c(-1, -2)]
# bbch_all_result$Country <- as.numeric(bbch_all_result$Country)
# bbch_all_result$Type <- as.numeric(bbch_all_result$Type)
# bbch_all_result$date <- as.numeric(bbch_all_result$date)

"STOP"

bbch_all_result[, -c(1)] <- lapply(bbch_all_result[, -c(1)], as.numeric)

bbch_all_result <- bbch_all_result[, c(-1, -2, -10)]
bbch_all_result <- bbch_all_result[, c(-10)]
bbch_all_result_edit <- bbch_all_result[, c(-1)]

M = cor(bbch_all_result_edit, method = "spearman")
testRes = cor.mtest(M, conf.level = 0.95)

p <- testRes$p

write.xlsx(M, "Plots/corrplots/non_normal/data_corr_Drought_allBBCHsep_sequenced_corrplot_lower_hclust_weight_seed.xlsx")
write.xlsx(p, "Plots/corrplots/non_normal/data_corr_Drought_allBBCHsep_sequenced_corrplot_lower_hclust_weight_seed_pvalues.xlsx")

write.xlsx(bbch_all_result, "Plots/corrplots/data_corr_phenotypes.xlsx")

pdf("Plots/corrplots/non_normal/Drought_allBBCHsep_sequenced_corrplot_lower_hclust_weight_seed_size.pdf", height=7, width = 7)
corrplot(M, p.mat = testRes$p, method = 'color', diag = T, type = 'lower',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20', order = 'hclust', outline = T, 
         col = COL2('RdYlBu'), tl.col = 'black') #Col1 is for sequential colours, Col2 is for diverging
dev.off()

library

### Correlation with groups not numbers -------------
#install.packages("rcompanion")

library(rcompanion)
#bbch_all_result$Admixed_K6_C

bbch_all_result$D <- bbch_all_result$GH
bbch_all_result$D <- gsub('Indeterminate_climbing', 'Indeterminate', bbch_all_result$D)
bbch_all_result$D <- gsub('Indeterminate_bush', 'Indeterminate', bbch_all_result$D)
bbch_all_result$D <- gsub('Determinate_bush', 'Determinate', bbch_all_result$D)

bbch_all_result$K2 <- bbch_all_result$Admixed_K6_C
bbch_all_result$K2 <- gsub('M1', 'Mesoamerican', bbch_all_result$K2)
bbch_all_result$K2 <- gsub('M2', 'Mesoamerican', bbch_all_result$K2)
bbch_all_result$K2 <- gsub('Admx_M', 'Mesoamerican', bbch_all_result$K2)
bbch_all_result$K2 <- gsub('A1', 'Andean', bbch_all_result$K2)
bbch_all_result$K2 <- gsub('C1', 'Andean', bbch_all_result$K2)
bbch_all_result$K2 <- gsub('C2', 'Andean', bbch_all_result$K2)
bbch_all_result$K2 <- gsub('C-EP', 'Andean', bbch_all_result$K2)
bbch_all_result$K2 <- gsub('Admx_A', 'Andean', bbch_all_result$K2)

colnames(bbch_all_result)[20] <- c("Major Seed Colour")
colnames(bbch_all_result)[21] <- c("SS")
colnames(bbch_all_result)[5] <- c("K6")
colnames(bbch_all_result)[24] <- c("PS")


discrete_df <- bbch_all_result[, c(1, 5:6, 20:26)] %>% unique()
discrete_df <- discrete_df[, -c(1)]


tbl <- table(discrete_df$Strategy, discrete_df$Admixed_K6_C)
tbl

cramerV(tbl)

chisq.test(tbl, simulate.p.value = TRUE, B = 10000)

df <- discrete_df %>%
  select(where(~ is.factor(.) || is.character(.))) %>%
  mutate(across(everything(), ~ factor(., exclude = NULL)))

var_names <- names(df)
n <- length(var_names)

cramer_matrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(var_names, var_names))
p_matrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(var_names, var_names))

for (i in 1:n) {
  for (j in i:n) {
    x <- df[[i]]
    y <- df[[j]]
    
    # Skip if not enough levels
    tbl <- table(x, y)
    if (any(dim(tbl) < 2)) next
    
    # Cramér’s V
    v <- tryCatch(cramerV(tbl), error = function(e) NA)
    
    # Chi-squared p-value
    p <- tryCatch(
      chisq.test(tbl, simulate.p.value = TRUE, B = 10000)$p.value,
      error = function(e) NA
    )
    
    # Fill both sides of the symmetric matrix
    cramer_matrix[i, j] <- v
    cramer_matrix[j, i] <- v
    p_matrix[i, j] <- p
    p_matrix[j, i] <- p
  }
  
  # Set diagonal
  cramer_matrix[i, i] <- 1
  p_matrix[i, i] <- NA
}

write.xlsx(cramer_matrix, "Plots/corrplots/discretevscon/chisq/cramer_matrix_discrete_seedsize.xlsx")
write.xlsx(p_matrix, "Plots/corrplots/discretevscon/chisq/cramer_pmatrix_discrete_seedsize.xlsx")

p_matrix[is.na(p_matrix)] <- 0

pdf("Plots/corrplots/discretevscon/chisq/Drought_allBBCHsep_discrete_rec.pdf", height=7, width = 7)
corrplot(cramer_matrix, p.mat = p_matrix, method = 'color', diag = T, type = 'lower',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20', order = 'hclust', outline = T, 
         col = COL2('RdYlBu'), col.lim = c(0, 1), tl.col = 'black') #Col1 is for sequential colours, Col2 is for diverging
dev.off()

#Discrete vs continuous, and continuous vs continuous 
library(rstatix)
#install.packages('flextable')
library('flextable')

colnames(bbch_all_result)[2] <- c("Date")
colnames(bbch_all_result)[14] <- "E apparent"
colnames(bbch_all_result)[15] <- "Fm'"
colnames(bbch_all_result)[13] <- "Gsw"
colnames(bbch_all_result)[12] <- "RH"

#which(colnames(bbch_all_result)=="RH")

df <- bbch_all_result[, c(-1, -21, -20, -7)] # For paper
df <- bbch_all_result[, c(-1, -4, -7)] #For thesis


#Last drought date and recovery date 
#bbch_all_result_30.08 <- bbch_all_result[grepl("date30.08", bbch_all_result$Date), ]
#bbch_all_result_23.08 <- bbch_all_result[grepl("date23.08", bbch_all_result$Date), ]

#df <- bbch_all_result_23.08[, -c(1,2)]
#df <- bbch_all_result_30.08[, -c(1,2)]

get_assoc_matrix_limited <- function(df) {
  df <- df %>%
    mutate(across(where(is.character), as.factor))  # Ensure consistent factors
  
  var_names <- names(df)
  n <- length(var_names)
  
  assoc_matrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(var_names, var_names))
  p_matrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(var_names, var_names))
  
  for (i in 1:n) {
    for (j in i:n) {
      x <- df[[i]]
      y <- df[[j]]
      name1 <- var_names[i]
      name2 <- var_names[j]
      
      is_x_num <- is.numeric(x)
      is_y_num <- is.numeric(y)
      is_x_fac <- is.factor(x)
      is_y_fac <- is.factor(y)
      
      assoc <- NA
      pval <- NA
      
      # Case 1: Continuous × Continuous (Spearman)
      if (is_x_num && is_y_num) {
        test <- suppressWarnings(cor.test(x, y, method = "spearman"))
        assoc <- test$estimate
        pval <- test$p.value
        
        # Case 2: Nominal × Continuous (Kruskal-Wallis + epsilon squared)
      } else if ((is_x_fac && is_y_num && nlevels(x) > 2) ||
                 (is_y_fac && is_x_num && nlevels(y) > 2)) {
        
        cont <- if (is_x_num) x else y
        disc <- if (is_x_fac) x else y
        df_sub <- data.frame(group = disc, value = cont)
        
        # Test + effect size
        test <- kruskal.test(value ~ group, data = df_sub)
        eff <- kruskal_effsize(df_sub, value ~ group)$effsize
        
        assoc <- eff
        pval <- test$p.value
        
        # Case 3: Discrete × Discrete → 0
      } else if (is_x_fac && is_y_fac) {
        assoc <- 0
        pval <- NA
      }
      
      # Fill both sides of symmetric matrices
      assoc_matrix[i, j] <- assoc
      assoc_matrix[j, i] <- assoc
      p_matrix[i, j] <- pval
      p_matrix[j, i] <- pval
    }
  }
  
  list(association = assoc_matrix, p_value = p_matrix)
}

result <- get_assoc_matrix_limited(df)
assoc_matrix <- result$association
p_matrix <- result$p_value
#p_matrix[is.na(p_matrix)] <- 0
assoc_matrix[is.na(assoc_matrix)] <- 0

pdf("Plots/corrplots/discretevscon/Drought_allBBCHsep_convdisc_new_thesis3.pdf", height=8, width = 8)
corrplot(assoc_matrix, p.mat = p_matrix, method = 'color', diag = T, type = 'lower',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20', order = 'hclust', outline = T, 
         col = COL2('RdYlBu'), tl.col = 'black') #Col1 is for sequential colours, Col2 is for diverging

dev.off()

write.xlsx(assoc_matrix, "Plots/corrplots/discretevscon/Drought_allBBCHsep_convdisc.xlsx")
write.xlsx(p_matrix, "Plots/corrplots/discretevscon/Drought_allBBCHsep_convdisc_pmatrirx.xlsx")


### Separate by date 
bbch_all_result_10.08 <- bbch_all_result[grepl("date10.08", bbch_all_result$Date), ]
bbch_all_result_18.08 <- bbch_all_result[grepl("date18.08", bbch_all_result$Date), ]
bbch_all_result_23.08 <- bbch_all_result[grepl("date23.08", bbch_all_result$Date), ]
bbch_all_result_30.08 <- bbch_all_result[grepl("date30.08", bbch_all_result$Date), ]
bbch_all_result_28.07 <- bbch_all_result[grepl("date28.07", bbch_all_result$Date), ]

bbch_all_result_10.08 <- bbch_all_result_10.08[, c(-1, -2)]
bbch_all_result_18.08 <- bbch_all_result_18.08[, c(-1, -2)]
bbch_all_result_23.08 <- bbch_all_result_23.08[, c(-1, -2)]
bbch_all_result_30.08 <- bbch_all_result_30.08[, c(-1, -2)]
bbch_all_result_28.07 <- bbch_all_result_28.07[, c(-1, -2)]


M = cor(bbch_all_result_23.08, method = "spearman")
testRes = cor.mtest(M, conf.level = 0.95)

pdf("Plots/corrplots/non_normal/Drought_23.08_sequenced_BBCHsep_corrplot_lower_hclust_weight_country_strat.pdf", height=7, width = 7)
corrplot(M, p.mat = testRes$p, method = 'color', diag = T, type = 'lower',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20', order = 'hclust', outline = T, 
         col = COL2('RdYlBu'), tl.col = 'black') #Col1 is for sequential colours, Col2 is for diverging
dev.off()

### Histograms for strategies ---------------------------------------
bbch_all_result_30.08$`Weight all foliar+pods`
bbch_all_result$Date <- gsub('date10.08', 'W2 Drought', bbch_all_result$Date)
bbch_all_result$Date <- gsub('date30.08', 'W5 Recovery', bbch_all_result$Date)
bbch_all_result$Date <- gsub('date28.07', 'W1 Drought', bbch_all_result$Date)
bbch_all_result$Date <- gsub('date23.08', 'W4 Drought', bbch_all_result$Date)
bbch_all_result$Date <- gsub('date18.08', 'W3 Drought', bbch_all_result$Date)
#bbch_all_result$Date <- as.Date(bbch_all_result$Date, format = "%d/%m/%Y")

pdf("Plots/Hist_pod_weight.pdf", height=8, width = 8)
p <- ggplot(bbch_all_result_30.08, aes(x=`Weight pods`)) + geom_histogram(bins = 15, color = "black") + 
  facet_grid(~Strategy) + xlab("Pod weight (g)") + ylab("Count") + theme_bw() + 
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12))
p
dev.off()

p1 <- ggplot(bbch_all_result_30.08, aes(x=`Weight all foliar+pods`)) + geom_histogram(bins = 15, color = "black") + 
  facet_grid(~Strategy) + xlab("Foliar weight (g)") + ylab("Count") + theme_bw() + 
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle = 0))
p1

p2 <- ggplot(bbch_all_result, aes(x=Date, y = BBCH, group = Date)) + geom_violin() +
  facet_grid(~Strategy) + xlab("Week") + ylab("BBCH") + theme_bw() +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))
p2

p3 <- ggplot(bbch_all_result, aes(x=Date, y = ETR, group = Date)) + geom_violin() +
  facet_grid(~Strategy) + xlab("Week") + theme_bw() +
  ylab(expression(ETR ~ (mu * mol ~ m^{-2} ~ s^{-1}))) +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))
p3

library(ggpubr)

pdf("Plots/Hist_violin_test.pdf", height=10, width = 12)
pall <- ggarrange(p, p1, p2, p3, nrow = 4, ncol = 1, labels = "AUTO", 
                  heights = c(1, 1, 1.4, 1.4))
pall
dev.off()

## For recovery 
# Tleaf, ETR, PhiPS2, E_apparent

Tleaf <- bbch_all_result[, c(1, 2, 7, 10)]
Tleaf_wide <- reshape(Tleaf, idvar = c("Accession", "Recovery"), timevar = "Date", direction = "wide")
Tleaf_wide$Diff_W4_W5 <- Tleaf_wide$`Tleaf.W4 Drought` -  Tleaf_wide$`Tleaf.W5 Recovery`

p <- ggplot(Tleaf_wide, aes(x=Recovery, y = Diff_W4_W5, group = Recovery)) + geom_violin() +
  geom_jitter(alpha = 0.5, position = position_jitter(width = 0.3)) +
   xlab("Recovery classification") + ylab("Difference in Tleaf (°C) \nbetween W4 and W5 ") + theme_bw() +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) +
  theme(axis.text.x = element_text())
p



### Try violin plots for Recovery ------------
# Tleaf, ETR, PhiPS2, E_apparent

bbch_all_result_lasttwo <- bbch_all_result[
  bbch_all_result$Date %in% c("W4 Drought", "W5 Recovery"),
]

write.xlsx(bbch_all_result_lasttwo, "all_data_BBCH_results_lasttwodates.xlsx")

bbch_all_result_last <- bbch_all_result[
  bbch_all_result$Date %in% c("W5 Recovery"),
]

bbch_all_result_fourth <- bbch_all_result[
  bbch_all_result$Date %in% c("W4 Drought"),
]

median(bbch_all_result_fourth$`Fm'`)
median(bbch_all_result_last$`Fm'`)

y_pos <- max(bbch_all_result_lasttwo$ETR, na.rm = TRUE) * 1.05

p <- ggplot(bbch_all_result_lasttwo, aes(x=Date, y = ETR, group = Date)) + 
  geom_violin() + 
  geom_jitter(alpha = 0.5, position = position_jitter(width = 0.2)) +
  xlab("Week") + ylab(expression(ETR ~ (mu * mol ~ m^{-2} ~ s^{-1}))) +
  theme_bw() +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) +
  theme(axis.text.x = element_text()) +
  
  ggsignif::geom_signif(
    comparisons = list(c("W4 Drought", "W5 Recovery")),
    test = "wilcox.test",
    map_signif_level = TRUE
  )

p


p1 <- ggplot(bbch_all_result_lasttwo, aes(x=Date, y = Tleaf, group = Date)) + 
  geom_violin() + 
  geom_jitter(alpha = 0.5, position = position_jitter(width = 0.2)) +
  xlab("Week") + ylab("Tleaf (°C)") +
  theme_bw() +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) +
  theme(axis.text.x = element_text()) +
  
  ggsignif::geom_signif(
    comparisons = list(c("W4 Drought", "W5 Recovery")),
    test = "wilcox.test",
    map_signif_level = TRUE)
p1

#bbch_all_result_lasttwo$`E apparent`

p2 <- ggplot(bbch_all_result_lasttwo, aes(x=Date, y = `E apparent`, group = Date)) + 
  geom_violin() + 
  geom_jitter(alpha = 0.5, position = position_jitter(width = 0.2)) +
  xlab("Week") + ylab(expression(italic(E) * " apparent" ~ (mmol ~ m^{-2} ~ s^{-1}))) +
  theme_bw() +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) +
  theme(axis.text.x = element_text()) +
  
  ggsignif::geom_signif(
    comparisons = list(c("W4 Drought", "W5 Recovery")),
    test = "wilcox.test",
    map_signif_level = TRUE)
p2

p3 <- ggplot(bbch_all_result_lasttwo, aes(x=Date, y = PhiPS2, group = Date)) + 
  geom_violin() + 
  geom_jitter(alpha = 0.5, position = position_jitter(width = 0.2)) +
  xlab("Week") + ylab(expression(PhiPS2 ~ (mu * mol ~ m^{-2} ~ s^{-1}))) +
  theme_bw() +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) +
  theme(axis.text.x = element_text()) +
  
  ggsignif::geom_signif(
    comparisons = list(c("W4 Drought", "W5 Recovery")),
    test = "wilcox.test",
    map_signif_level = TRUE)
p3

p4 <- ggplot(bbch_all_result_lasttwo, aes(x=Date, y = BBCH, group = Date)) + 
  geom_violin() + 
  geom_jitter(alpha = 0.5, position = position_jitter(width = 0.2)) +
  xlab("Week") + ylab("BBCH") +
  theme_bw() +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) +
  theme(axis.text.x = element_text()) +
  ggsignif::geom_signif(
    comparisons = list(c("W4 Drought", "W5 Recovery")),
    test = "wilcox.test",
    map_signif_level = TRUE)
p4

p5 <- ggplot(bbch_all_result_lasttwo, aes(x=Date, y = Fs, group = Date)) + 
  geom_violin() + 
  geom_jitter(alpha = 0.5, position = position_jitter(width = 0.2)) +
  xlab("Week") + ylab("Fs") +
  theme_bw() +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) +
  theme(axis.text.x = element_text()) +
  ggsignif::geom_signif(
    comparisons = list(c("W4 Drought", "W5 Recovery")),
    test = "wilcox.test",
    map_signif_level = TRUE)
p5

p6 <- ggplot(bbch_all_result_lasttwo, aes(x=Date, y = RH, group = Date)) + 
  geom_violin() + 
  geom_jitter(alpha = 0.5, position = position_jitter(width = 0.2)) +
  xlab("Week") + ylab("RH (%)") +
  theme_bw() +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) +
  theme(axis.text.x = element_text()) +
  ggsignif::geom_signif(
    comparisons = list(c("W4 Drought", "W5 Recovery")),
    test = "wilcox.test",
    map_signif_level = TRUE)
p6

p7 <- ggplot(bbch_all_result_lasttwo, aes(x=Date, y = Gsw, group = Date)) + 
  geom_violin() + 
  geom_jitter(alpha = 0.5, position = position_jitter(width = 0.2)) +
  xlab("Week") + ylab(expression(Gsw ~ (mol ~ m^{-2} ~ s^{-1}))) + 
  theme_bw() +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) +
  theme(axis.text.x = element_text()) +
  ggsignif::geom_signif(
    comparisons = list(c("W4 Drought", "W5 Recovery")),
    test = "wilcox.test",
    map_signif_level = TRUE)
p7

bbch_all_result_lasttwo$`Fm'`

p8 <- ggplot(bbch_all_result_lasttwo, aes(x=Date, y = `Fm'`, group = Date)) + 
  geom_violin() + 
  geom_jitter(alpha = 0.5, position = position_jitter(width = 0.2)) +
  xlab("Week") + ylab("Fm'") + 
  theme_bw() +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) +
  theme(axis.text.x = element_text()) +
  ggsignif::geom_signif(
    comparisons = list(c("W4 Drought", "W5 Recovery")),
    test = "wilcox.test",
    map_signif_level = TRUE)
p8

p9 <- ggplot(bbch_all_result_lasttwo, aes(x=Date, y = VPDleaf, group = Date)) + 
  geom_violin() + 
  geom_jitter(alpha = 0.5, position = position_jitter(width = 0.2)) +
  xlab("Week") + ylab("VPDleaf (kPa)") + 
  theme_bw() +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) +
  theme(axis.text.x = element_text()) +
  ggsignif::geom_signif(
    comparisons = list(c("W4 Drought", "W5 Recovery")),
    test = "wilcox.test",
    map_signif_level = TRUE)
p9

pdf("Plots/Recovery_violin_test.pdf", height=15, width = 12)
p_all <- ggarrange(p, p1, p2, p3, p4, p5, p6, p7, p8, 
                   ncol = 3, nrow = 3, 
                   labels = "AUTO")
p_all
dev.off()

### Histograms for drought related traits ----------------
#Don't remove controls
harvest <- read_excel("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Bean project/23.03-Field-trial-NIAB/KDJ Bean panel randomisation.xlsx", sheet = "harvest_edit")
all_df <- read.xlsx("Plots/corrplots/all_dates_licor.xlsx")

K6_clustering <- read_excel("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Bean project/STATS/22.01 Variant calling - GATK/22.11-ADMIXTURE-LD0.5/Q_files/K6_Q_files/K6_clustering_PCA-phenotypes over 0.7 ADMX.xlsx")
photosens <- K6_clustering[, c(2, 9, 10, 37, 47)]
colnames(photosens)[1] <- c("Accession")
colnames(photosens)[5] <- c("photo_sens")

harvest <- merge(harvest, photosens, by = "Accession") %>%
  unique()

harvest$Accession_T <- paste(harvest$Accession, harvest$Treat, sep = "_") 
harvest <- harvest[, c(23, 6, 10, 12:14, 22, 24:26)]

harvest[, -c(1,2, 3, 7:10)] <- lapply(harvest[, -c(1,2, 3, 7:10)], as.numeric)
harvest <- aggregate(. ~ Accession_T + GH + Treat + photo_sens + Admixed_A_K2 + Type + Country, 
                     data = harvest, 
                     FUN = mean, 
                     na.rm = TRUE)


write.xlsx(harvest, "mean_weights_pods_foliar_all_correct.xlsx")

all_join <- join(all_df, harvest, by = "Accession_T") %>%
  unique()

all_join[, -c(1, 11:17)] <- lapply(all_join[, -c(1, 11:17)], as.numeric)

result <- aggregate(. ~ Accession_T + GH + date + Treat + photo_sens + Admixed_A_K2 + Type + Country, 
                    data = all_join, 
                    FUN = mean, 
                    na.rm = TRUE)

accession_list <- c("G50516A_D", "G50516O_D", "G50516D_D", "G51285_D", "G51285H_D", "G8152_D", "G916_D", "G50997_D",
                    "G50516A_C", "G50516O_C", "G50516D_C", "G51285_C", "G51285H_C", "G8152_C", "G916_C", "G50997_C"
                    )  
result <- result[!result$Accession_T %in% accession_list, ]

#With BBCH
date30.08 <- read.xlsx("30.08.23.merged.xlsx")
date23.08 <- read.xlsx("23.08.23.merged.xlsx")
date18.08 <- read.xlsx("18.08.23.merged.xlsx")
date10.08 <- read.xlsx("10.08.23.merged.xlsx")
date28.07 <- read.xlsx("28.07.23.merged.xlsx")

date30.08$Accession_T <- paste(date30.08$Accession, date30.08$Treat, sep = "_")
date23.08$Accession_T <- paste(date23.08$Accession, date23.08$Treat, sep = "_")
date18.08$Accession_T <- paste(date18.08$Accession, date18.08$Treat, sep = "_")
date10.08$Accession_T <- paste(date10.08$Accession, date10.08$Treat, sep = "_")
date28.07$Accession_T <- paste(date28.07$Accession, date28.07$Treat, sep = "_")

date30.08 <- date30.08[, c(111, 106)]
date23.08 <- date23.08[, c(111, 106)]
date18.08 <- date18.08[, c(111, 106)]
date10.08 <- date10.08[, c(111, 106)]
date28.07 <- date28.07[, c(111, 106)]

date30.08$date <- c("date30.08")
date23.08$date <- c("date23.08")
date18.08$date <- c("date18.08")
date10.08$date <- c("date10.08")
date28.07$date <- c("date28.07")

colnames(date30.08)[2] <- c("BBCH")
colnames(date23.08)[2] <- c("BBCH")
colnames(date18.08)[2] <- c("BBCH")
colnames(date10.08)[2] <- c("BBCH")
colnames(date28.07)[2] <- c("BBCH")

df_list <- list(date30.08, date23.08, date18.08, date10.08, date28.07)

BBCH <- bind_rows(df_list)

BBCH_res <- aggregate(. ~ Accession_T + date, 
                      data = BBCH, 
                      FUN = mean, 
                      na.rm = TRUE)

BBCH_res <- BBCH_res[!BBCH_res$Accession_T %in% accession_list, ]

bbch_all_result <- join(result, BBCH_res, by = c("Accession_T", "date")) %>%
  unique()
bbch_all_result$Treat_GH <- paste(bbch_all_result$GH, bbch_all_result$Treat, sep = "_")
colnames(bbch_all_result)[2] <- c("Date")
bbch_all_result$Date <- gsub('date28.07', '28/07/2023', bbch_all_result$Date)
bbch_all_result$Date <- gsub('date10.08', '10/08/2023', bbch_all_result$Date)
bbch_all_result$Date <- gsub('date18.08', '18/08/2023', bbch_all_result$Date)
bbch_all_result$Date <- gsub('date23.08', '23/08/2023', bbch_all_result$Date)
bbch_all_result$Date <- gsub('date30.08', '30/08/2023', bbch_all_result$Date)
level_order <- c("28/07/2023", "10/08/2023", "18/08/2023", "23/08/2023", "30/08/2023")
bbch_all_result$Date <- factor(bbch_all_result$Date, levels=level_order)

write.xlsx(bbch_all_result, "mean_weights_licor_bbch_date_long_country_pop.xlsx")
bbch_all_result <- read.xlsx("mean_weights_licor_bbch_date_long_country_pop.xlsx")

bbch_all_result_drought <- bbch_all_result[!grepl("_C", bbch_all_result$Accession_T), ]
bbch_all_result_drought_last <- bbch_all_result_drought[grepl("30/08/2023", bbch_all_result_drought$Date), ]


facet_labels <- c(
  "Determinate_bush_C" = "Control DB",
  "Determinate_bush_D" = "Drought DB",
  "Indeterminate_bush_D" = "Drought IB",
  "Indeterminate_climbing_D" = "Drought IC"
)

pdf("Plots/corrplots/phenos_by_date/All_drought_hist_PhiPS2_GH.pdf", height=6, width = 14)
p1 <- ggplot(bbch_all_result, aes(x=PhiPS2, fill = Date)) + geom_histogram(bins = 10) + 
  facet_grid(~Treat_GH, labeller = labeller(Treat_GH = as_labeller(facet_labels))) + 
  xlab("gsw (mol H2O m-2 s-1)") + ylab("Count") + theme_bw() +
  scale_fill_brewer(palette = "Dark2") +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12))
p1
dev.off()

##Quantile values
prop <- read.xlsx("mean_weights_pods_foliar_all_correct.xlsx")
prop_D <- prop[grepl("D", prop$Treat),]
quantile(prop_D$Proportion, na.rm = T, probs = c(0.25, 0.5, 0.75, 0.95))

options(scipen = 999)
harvesting <- bbch_all_result[, c(1, 4, 18,19, 20)] %>%
  unique()
harvesting <- harvesting[grepl("D", harvesting$Treat),]

quantile(bbch_all_result_drought_last$PhiPS2, probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))

E <- bbch_all_result_drought_last[, c(1:4, 15, 31)]

bbch_all_result_drought <- bbch_all_re
Weight.podsbbch_all_result_drought <- bbch_all_result[grepl("D", bbch_all_result$Treat),]

bbch_all_result_last_drought <- bbch_all_result_drought[grepl("30/08/2023", bbch_all_result_drought$Date), ]
bbch_all_result_first_drought <- bbch_all_result_drought[grepl("28/07/2023", bbch_all_result_drought$Date), ]
bbch_all_result_third_drought <- bbch_all_result_drought[grepl("18/08/2023", bbch_all_result_drought$Date), ]
bbch_all_result_fourth_drought <- bbch_all_result_drought[grepl("23/08/2023", bbch_all_result_drought$Date), ]

quantile(bbch_all_result_first_drought$gsw)

#Difference between last two dates for recovery information
Tleaf <- bbch_all_result_drought[, c(1:3,11)]
Tleaf_wide <- reshape(Tleaf, idvar = c("Accession_T", "GH"), timevar = "Date", direction = "wide")
Tleaf_wide$Diff_fourth_fifth <- Tleaf_wide$`Tleaf.23/08/2023` - Tleaf_wide$`Tleaf.30/08/2023`

quantile(Tleaf_wide$Diff_fourth_fifth, probs = 0.95)

ETR <- bbch_all_result_drought[, c(1:3,9)]
ETR_wide <- reshape(ETR, idvar = c("Accession_T", "GH"), timevar = "Date", direction = "wide")
ETR_wide$Diff_fourth_fifth <- ETR_wide$`ETR.23/08/2023` - ETR_wide$`ETR.30/08/2023`

quantile(ETR_wide$Diff_fourth_fifth, probs = 0.1)

PhiPS2 <- bbch_all_result_drought[, c(1:3,17)]
PhiPS2_wide <- reshape(PhiPS2, idvar = c("Accession_T", "GH"), timevar = "Date", direction = "wide")
PhiPS2_wide$Diff_fourth_fifth <- PhiPS2_wide$`PhiPS2.23/08/2023` - PhiPS2_wide$`PhiPS2.30/08/2023`

quantile(PhiPS2_wide$Diff_fourth_fifth, probs = 0.1)

E <- bbch_all_result_drought[, c(1:3,15)]
E_wide <- reshape(E, idvar = c("Accession_T", "GH"), timevar = "Date", direction = "wide")
E_wide$Diff_fourth_fifth <- E_wide$`E_apparent.23/08/2023` - E_wide$`E_apparent.30/08/2023`

quantile(E_wide$Diff_fourth_fifth, probs = 0.05)

### Violin plots --------------### ViolE_apparentin plots --------------

bbch_all_result_last <- bbch_all_result[grepl("30/08/2023", bbch_all_result$Date), ]

?ggbetweenstats
#non-parametric for non-normal
#For stricter control of false positives → Holm also default
pdf("Plots/corrplots/phenos_by_date/gsw_all_GH_violin.pdf", height=8, width = 8)
p <- ggstatsplot::ggbetweenstats(data = bbch_all_result_last,
                                     x = Treat_GH,
                                     y = gsw,
                                     type = "np",
                                     pairwise.display = "significant",
                                     p.adjust.method = "holm",
                                     results.subtitle = FALSE,
                                     centrality.point.args = list(size = 4, color = "darkred"),
                                     violin.args = list(width = 1, alpha = 0.2, na.rm = TRUE),
                                     ggsignif.args = list(textsize = 3, tip_length = 0.01, na.rm = TRUE),
                                     # mean.ci = TRUE,
                                     pairwise.comparisons = TRUE,
                                     boxplot.args = list(width = 0),
                                     bf.message = FALSE,
                                     xlab = "Treatment and Growth Habit", 
                                     ylab = "gsw (mmol m−2 s−1 )", #ETR (µmol s-1), Tleaf (°C), gsw (mmol m−2 s−1 )
                                     digits = "signif",
                                     ggtheme = ggplot2::theme_bw()
) +
  scale_x_discrete(labels = c(
    "Determinate_bush_C" = "Control \nDeterminate Bush", 
    "Determinate_bush_D" = "Drought \nDeterminate Bush", 
    "Indeterminate_bush_D" = "Drought \nIndeterminate Bush", 
    "Indeterminate_climbing_D" = "Drought \nIndeterminate Climbing"
  )) #+
# theme(
#   axis.text.x = element_text(angle = 45, hjust = 1))
p
dev.off()

bbch_all_result_drought_lasttwo$PhiPS2

pdf("Plots/corrplots/phenos_by_date/Tleaf_date_violin.pdf", height=8, width = 8)
p3 <- ggstatsplot::ggbetweenstats(data = bbch_all_result_drought_lasttwo,
                                 x = Date,
                                 y = Tleaf,
                                 type = "np",
                                 pairwise.display = "significant",
                                 p.adjust.method = "holm",
                                 results.subtitle = FALSE,
                                 centrality.point.args = list(size = 3, color = "darkred"),
                                 violin.args = list(width = 1, alpha = 0.2, na.rm = TRUE),
                                 ggsignif.args = list(textsize = 3, tip_length = 0.01, na.rm = TRUE, 
                                                      map_signif_level = TRUE),
                                 # mean.ci = TRUE,
                                 pairwise.comparisons = TRUE,
                                 boxplot.args = list(width = 0.2),
                                 bf.message = FALSE,
                                 xlab = "Date", 
                                 ylab = "Tleaf (°C)", #ETR (µmol s-1), Tleaf (°C), gsw (mmol m−2 s−1 ), ylab(expression(PhiPS2 ~ (mu * mol ~ m^{-2} ~ s^{-1}))), expression(E_apparent ~ (mol ~ m^{-2} ~ s^{-1}))
                                 digits = "signif",
                                 ggtheme = ggplot2::theme_bw(), 
                                 label.args = list(n.text = FALSE) 
)  +
  scale_color_manual(values = rep("black", length(unique(bbch_all_result_drought_lasttwo$Date)))) +
  scale_fill_manual(values = rep("grey80", length(unique(bbch_all_result_drought_lasttwo$Date))))
# p +
#  ggplot2::geom_point(aes(color = Strategy), position = position_jitter(width = 0.1), size = 2, data = bbch_all_result_drought_lasttwo) +
#  ggplot2::scale_color_brewer(palette = "Dark2")

p3
dev.off()

### Scatter plots --------
library(ggrepel)
library(ggforce)

#bbch_all_result_drought$`Weight pods`

bbch_all_result <- read.xlsx("mean_weights_licor_bbch_date_long_country_pop.xlsx")
bbch_all_result[1, 1] <- "BONNEBOUCHE_D"

bbch_all_result_2308 <- bbch_all_result[grepl("23/08/2023", bbch_all_result$Date), ]

bbch_all_result_drought <- bbch_all_result[!grepl("_C", bbch_all_result$Accession_T), ]
bbch_all_result_drought_last <- bbch_all_result_drought[grepl("30/08/2023", bbch_all_result_drought$Date), ]
bbch_all_result_drought_2308 <- bbch_all_result_drought[grepl("23/08/2023", bbch_all_result_drought$Date), ]

bbch_all_result_drought_lasttwo <- bbch_all_result_drought[
  grepl("23/08/2023|30/08/2023", bbch_all_result_drought$Date),
]

bbch_all_result_drought <- bbch_all_result_drought[, -c(22:23)]
bbch_all_result_drought_last <- bbch_all_result_drought_last[, -c(22:23)]
bbch_all_result_drought_2308 <- bbch_all_result_drought_2308[, -c(22:23)]
bbch_all_result_drought_lasttwo <- bbch_all_result_drought_lasttwo[, -c(22:23)]

# bbch_all_result_drought_last$Accession <- bbch_all_result_drought_last$Accession_T
# bbch_all_result_drought_last <- bbch_all_result_drought_last %>% separate(Accession, c("Accession", NA))
# 
# bbch_all_result_drought$Accession <- bbch_all_result_drought$Accession_T
# bbch_all_result_drought <- bbch_all_result_drought %>% separate(Accession, c("Accession", NA))
# 
# bbch_all_result$Accession <- bbch_all_result$Accession_T
# bbch_all_result <- bbch_all_result %>% separate(Accession, c("Accession", "Treat"))

# bbch_all_result_drought <- bbch_all_result_drought %>%
#   mutate(filtered_accession = ifelse((`Weight all foliar+pods` > 140 | `Weight all foliar+pods` < 52) | `Weight pods` > 41, Accession, NA))

bbch_all_result_drought$GH <- gsub('Determinate_bush', 'Determinate bush', bbch_all_result_drought$GH)
bbch_all_result_drought$GH <- gsub('Indeterminate_bush', 'Indeterminate bush', bbch_all_result_drought$GH)
bbch_all_result_drought$GH <- gsub('Indeterminate_climbing', 'Indeterminate climbing', bbch_all_result_drought$GH)

# accession_list <- c("G12712", "G23462", "G50516S", "MADEIRAMARRON", "PHA13866", "PHA13872", "PHA14071")
# descriptor_list <- c("Stay-green", "Drought sensitive", "Drought sensitive", "Drought tolerant", "Drought tolerant", "Stay-green", "Drought tolerant")
# bbch_all_result_drought <- bbch_all_result_drought %>%
#   mutate(Descriptor = descriptor_list[match(Accession, accession_list)])

strategies <- read.xlsx("mean_weights_pods_foliar_all_correct.xlsx")
strategies <- strategies[, c(1, 8:17, 19)]

bbch_all_result_drought <- join(bbch_all_result_drought, strategies, by = "Accession_T") %>%
  unique()

bbch_all_result_drought_last <- join(bbch_all_result_drought_last, strategies, by = "Accession_T") %>%
  unique()
colnames(bbch_all_result_drought_last)[31] <- c("Strategy")

bbch_all_result_drought_2308 <- join(bbch_all_result_drought_2308, strategies, by = "Accession_T") %>%
  unique()
colnames(bbch_all_result_drought_2308)[31] <- c("Strategy")

bbch_all_result_drought_lasttwo <- join(bbch_all_result_drought_lasttwo, strategies, by = "Accession_T") %>%
  unique()
colnames(bbch_all_result_drought_lasttwo)[31] <- c("Strategy")

#bbch_all_result_drought_last$Weight.pods

str(bbch_all_result_drought_last)

pdf("Plots/scatterplots/All_drought_last_weightpods_weightfoliar_colstrat.pdf", height=5, width = 7)
p1 <- ggplot(bbch_all_result_drought_last, aes(x=`Weight.all.foliar+pods`, y = `Weight.pods`, colour = Strategy, fill = Strategy)) + 
  geom_point() + 
  xlab("Foliar weight (g)") + ylab("Pod weight (g)") + 
  theme_bw() +
  scale_fill_brewer(palette = "Dark2") +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) 
  # ggforce::geom_mark_ellipse(data = subset(bbch_all_result_drought_last, !(Strategy %in% c("No classification"))),
  #   aes(), alpha = 0, show.legend = FALSE, tol = 0.05,
  #   expand = unit(1, "mm"))
#+
  # geom_label_repel(aes(label = Accession, colour = Descriptor), size = 3,
  #                 box.padding = 0.8, #0.3 and 0.5 works well but not lines, 0.7 and 0.4 have lines
  #                 point.padding = 0.4, #0.4 is closest
  #                 nudge_y = 1, 
  #                 ylim = c(0,65)
  #                   ) +
  # scale_colour_discrete(na.translate = FALSE) 
p1
dev.off()

?geom_mark_ellipse

# bbch_all_result <- bbch_all_result %>%
#   mutate(Descriptor = descriptor_list[match(Accession, accession_list)])
# bbch_all_result <- bbch_all_result %>%
#   mutate(filtered_accession = accession_list[match(Accession, accession_list)])
# 
# bbch_all_result_last <- bbch_all_result[grepl("30/08/2023", bbch_all_result$Date), ]
# 
# bbch_all_result_drought_last <- bbch_all_result_drought_last %>%
#   mutate(Descriptor = descriptor_list[match(Accession, accession_list)])
# bbch_all_result_drought_last <- bbch_all_result_drought_last %>%
#   mutate(filtered_accession = accession_list[match(Accession, accession_list)])


pdf("Plots/scatterplots/All_drought_weightfoliar_gsw_2308_colgrouping.pdf", height=6, width = 8)
p3 <- ggplot(bbch_all_result_drought_2308, aes(x= ETR, y = PhiPS2, colour = Strategy)) + geom_point() + 
  xlab(expression(ETR ~ (mu * mol ~ m^{-2} ~ s^{-1}))) + ylab(expression(PhiPS2 ~ (mu * mol ~ m^{-2} ~ s^{-1}))) + 
  theme_bw() +
  scale_fill_brewer(palette = "Dark2") +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) +
  ggforce::geom_mark_ellipse(data = subset(bbch_all_result_drought_2308, !(Strategy %in% c("No classification", "Saver", "Prioritised yield"))),
                             aes(), alpha = 0.1, show.legend = FALSE, tol = 0.05,
                             expand = unit(1, "mm"))

#+
  # geom_label_repel(aes(label = Accession, colour = Descriptor), size = 3,
  #                  box.padding = 0.8, #0.3 and 0.5 works well but not lines, 0.7 and 0.4 have lines
  #                  point.padding = 0.4 #0.4 is closest
  # ) +
  # scale_colour_discrete(na.translate = FALSE) 
p3

dev.off()

pdf("Plots/scatterplots/All_drought_gsw_Tleaf_weightlabels_2308.pdf", height=6, width = 8)
p1 <- ggplot(bbch_all_result_drought_2308, aes(x=gsw, y = Tleaf, colour = Second.grouping)) + geom_point() + 
  xlab("gsw (mmol m−2 s−1)") + ylab("Tleaf (°C)") + 
  theme_bw() +
  scale_fill_brewer(palette = "Dark2") +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) +
  geom_label_repel(aes(label = Accession, colour = Descriptor), size = 3,
                  box.padding = 0.8, 
                  point.padding = 0.5
                  
  ) +
  scale_colour_discrete(na.translate = FALSE)
p1
dev.off()

#ETR (µmol s-1), Tleaf (°C), gsw (mmol m−2 s−1), Fm (\mu molm^{-2}s^{-1}\)


# bbch_all_result_drought <- bbch_all_result_drought %>%
#   mutate(filtered_accession = ifelse((Tleaf > 23 | Tleaf < 19) | gsw > 0.7, Accession, NA))

Tleaf_wide$`Tleaf.23/08/2023`
Tleaf_wide <- join(Tleaf_wide, strategies, by = "Accession_T")
colnames(Tleaf_wide)[18] <- c("Strategy")

pdf("Plots/scatterplots/All_drought_gsw_Tleaf.pdf", height=6, width = 8)
p4 <- ggplot(Tleaf_wide, aes(x=`Tleaf.23/08/2023`, y = `Tleaf.30/08/2023`, colour = Strategy)) + geom_point() + 
  xlab("Tleaf (°C) 23/08/2023") + ylab("Tleaf (°C) 30/08/2023") + 
  theme_bw() +
  scale_fill_brewer(palette = "Dark2") +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) +
 ggforce::geom_mark_ellipse(data = subset(Tleaf_wide,  !(Strategy %in% c("No classification", "DB control growth", "Saver", "Prioritised yield"))),
                            aes(), alpha = 0.1, show.legend = FALSE, tol = 0.05, 
                            expand = unit(1, "mm"))

#+
  # geom_label_repel(aes(label = filtered_accession, Descriptor), size = 3,
  #                  box.padding = 0.8, 
  #                  point.padding = 0.5
  #                  
  # ) +
  # scale_colour_discrete(na.translate = FALSE)
p4
dev.off()

#Fluorometer
# bbch_all_result_drought <- bbch_all_result_drought %>%
#   mutate(filtered_accession = ifelse((PhiPS2 > 0.63 | PhiPS2 < 0.2) | (ETR < 40 | ETR > 150), Accession, NA))

colnames(bbch_all_result_drought_2308)[31] <- c("Strategy")

bbch_all_result_drought_2308$gswsq <- bbch_all_result_drought_2308$gsw^2

#ETR and PhiPS2
#ylab(expression(PhiPS2 ~ (mu * mol ~ m^{-2} ~ s^{-1})))
#ylab(expression(gsw ~ (mmol ~ m^{-2} ~ s^{-1})))
#xlab(expression(ETR ~ (mu * mol ~ m^{-2} ~ s^{-1})))
#xlab(expression(E apparent ~ (mol ~ m^{-2} ~ s^{-1})))
pdf("Plots/scatterplots/All_drought_gsw_E_apparent_labels_2308_colstrat.pdf", height=5, width = 7)
p4 <- ggplot(bbch_all_result_drought_2308, aes(x=gswsq, y = E_apparent, colour = Strategy)) + 
  geom_point() + 
  xlab(expression(gsw ~ (mmol ~ m^{-2} ~ s^{-1}))) + ylab(expression(E_apparent ~ (mol ~ m^{-2} ~ s^{-1}))) + 
  theme_bw() +
  scale_fill_brewer(palette = "Dark2") +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) +
  geom_jitter(width = 0.001, height = 0)
#+
  # geom_label_repel(aes(label = Accession, colour = Descriptor), size = 3, #, colour = Descriptor
  #                  box.padding = 0.8,
  #                  point.padding = 0.5
  # 
  # ) +
  # scale_colour_discrete(na.translate = FALSE)
p4
dev.off()

##BBCH

BBCH_res <- bbch_all_result_drought[, c(1, 3, 21)]
BBCH_res[1, 1] <- "BONNEBOUCHE_D"

BBCH_wide <- reshape(BBCH_res, idvar = "Accession_T", timevar = "Date", direction = "wide")
#BBCH_wide$Accession <- BBCH_wide$Accession_T
#BBCH_wide <- BBCH_wide %>% separate(Accession, c("Accession", "Treat"))

# accession_list <- c("G12712", "G23462", "G50516S", "MADEIRAMARRON", "PHA13866", "PHA13872", "PHA14071")
# descriptor_list <- c("Stay-green", "Drought sensitive", "Drought sensitive", "Drought tolerant", "Drought tolerant", "Stay-green", "Drought tolerant")
# BBCH_wide <- BBCH_wide %>%
#   mutate(Descriptor = descriptor_list[match(Accession, accession_list)])
# 
# BBCH_wide <- BBCH_wide %>%
#   mutate(filtered_accession = accession_list[match(Accession, accession_list)])

#BBCH_wide$BBCH.date23.08

BBCH_wide$diff_fourth_fifth <- BBCH_wide$`BBCH.23/08/2023` - BBCH_wide$`BBCH.30/08/2023`

BBCH_wide <- join(BBCH_wide, strategies, by = "Accession_T") %>%
  unique()

colnames(BBCH_wide)[16] <- c("Strategy")

pdf("Plots/scatterplots/All_drought_BBCH_23_30_colstrat.pdf", height=5, width = 7)
p2 <- ggplot(BBCH_wide, aes(x=`BBCH.23/08/2023`, y = `BBCH.30/08/2023`, colour = Strategy)) + geom_point() + 
  xlab("BBCH 23/08/2023") + ylab("BBCH 30/08/2023") + 
  theme_bw() +
  scale_fill_brewer(palette = "Dark2") +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) 
  # ggforce::geom_mark_ellipse(data = subset(BBCH_wide,  !(Strategy %in% c("No classification"))),
  #                            aes(), alpha = 0.1, show.legend = FALSE, tol = 0.05,
  #                            expand = unit(1, "mm"))

#+
  # geom_label_repel(aes(label = Accession, colour = Descriptor), size = 3, #, colour = Descriptor
  #                  box.padding = 0.8, 
  #                  point.padding = 0.5
  #                  
  # ) +
  # scale_colour_discrete(na.translate = FALSE) +
  # scale_y_continuous(limits = c(20, 90))
p2
dev.off()

#### Differences for recovery 
bbch_all_result_drought[1, 1] <- "BONNEBOUCHE_D"
colnames(bbch_all_result_drought)[31] <- c("Strategy")

Tleaf_E <- bbch_all_result_drought[, c(1, 3, 11, 15, 31)]
PhiPS2_ETR <- bbch_all_result_drought[, c(1, 3, 9, 17, 31)]

strategies <- read.xlsx("mean_weights_pods_foliar_all_correct.xlsx")
strategies <- strategies[, c(1, 19)]
colnames(strategies)[2] <- c("Recovery (Difference between 23/08/23 and 30/08/2023)")

Tleaf_E <- Tleaf_E %>% filter(Date %in% c("23/08/2023", "30/08/2023"))
Tleaf_E_wide <- reshape(Tleaf_E, idvar = c("Accession_T", "Strategy"), timevar = "Date", direction = "wide")
Tleaf_E_wide$diff_Tleaf <- Tleaf_E_wide$`Tleaf.23/08/2023` -  Tleaf_E_wide$`Tleaf.30/08/2023`
Tleaf_E_wide$diff_E <- Tleaf_E_wide$`E_apparent.23/08/2023` -  Tleaf_E_wide$`E_apparent.30/08/2023`
Tleaf_E_wide <- left_join(Tleaf_E_wide, strategies, by = "Accession_T")

PhiPS2_ETR <- PhiPS2_ETR %>% filter(Date %in% c("23/08/2023", "30/08/2023"))
PhiPS2_ETR_wide <- reshape(PhiPS2_ETR, idvar = c("Accession_T", "Strategy"), timevar = "Date", direction = "wide")
PhiPS2_ETR_wide$diff_PhiPS2 <- PhiPS2_ETR_wide$`PhiPS2.23/08/2023` - PhiPS2_ETR_wide$`PhiPS2.30/08/2023`
PhiPS2_ETR_wide$diff_ETR <- PhiPS2_ETR_wide$`ETR.23/08/2023` - PhiPS2_ETR_wide$`ETR.30/08/2023`
PhiPS2_ETR_wide <- left_join(PhiPS2_ETR_wide, strategies, by = "Accession_T")


#ylab(expression(PhiPS2 ~ (mu * mol ~ m^{-2} ~ s^{-1})))
#ylab(expression(gsw ~ (mmol ~ m^{-2} ~ s^{-1})))
#xlab(expression(ETR ~ (mu * mol ~ m^{-2} ~ s^{-1})))
#xlab(expression(E apparent ~ (mol ~ m^{-2} ~ s^{-1})))
#Tleaf °C

PhiPS2_ETR_wide$`Recovery (Difference between 23/08/23 and 30/08/2023)`

pdf("Plots/scatterplots/All_drought_PhiPS2_ETR_apparent_23_30_colstrat.pdf", height=5, width = 7)
p6 <- ggplot(Tleaf_E_wide, aes(x=diff_E, y = diff_Tleaf, colour = `Recovery (Difference between 23/08/23 and 30/08/2023)`)) + 
  geom_point() + 
  xlab(expression(E_apparent ~ (mol ~ m^{-2} ~ s^{-1}))) + ylab("Tleaf °C") + 
  theme_bw() +
  scale_fill_brewer(palette = "Dark2") +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) +
  geom_jitter(width = 0.001, height = 0) +
  labs(colour = "Recovery (Difference between 23/08/23 and 30/08/2023)") 
  # ggforce::geom_mark_ellipse(data = subset(Tleaf_E_wide,  !(Strategy %in% c("No classification"))),
  #                            aes(), alpha = 0.1, show.legend = FALSE, tol = 0.05,
  #                            expand = unit(1, "mm"))

p6
dev.off()

p7 <- ggplot(PhiPS2_ETR_wide, aes(x=diff_ETR, y = diff_PhiPS2, colour = `Recovery (Difference between 23/08/23 and 30/08/2023)`)) + 
  geom_point() + 
  xlab(expression(ETR ~ (mu * mol ~ m^{-2} ~ s^{-1}))) + ylab(expression(PhiPS2 ~ (mu * mol ~ m^{-2} ~ s^{-1}))) + 
  theme_bw() +
  scale_fill_brewer(palette = "Dark2") +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) +
  geom_jitter(width = 0.001, height = 0) +
  labs(colour = "Recovery (Difference between 23/08/23 and 30/08/2023)") 
  # ggforce::geom_mark_ellipse(data = subset(PhiPS2_ETR_wide,  !(Strategy %in% c("No classification"))),
  #                            aes(), alpha = 0.1, show.legend = FALSE, tol = 0.05,
  #                            expand = unit(1, "mm"))
p7

### Put all plots together 

p_rec <- ggarrange(p6, p7, nrow = 1, ncol = 2, common.legend = T, legend = "bottom", 
                   labels = c("C", "D"), vjust = 1)
p_rec
?ggarrange

# p_other <- ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, common.legend = T, legend = "right", 
#                     labels = c("A", "B", "C", "D"), vjust = 1)

p_scat <- ggarrange(p1, p2, nrow = 1, ncol = 2, common.legend = T, legend = "bottom", 
                     labels = c("A", "B"), vjust = 1)
p_scat

pdf("Plots/scatterplots/All_plots_drought_scatters_bottomleg.pdf", height=8, width = 8)
p_all <- ggarrange(p_scat, p_rec, ncol = 1, nrow = 2, widths = c(1, 1), 
                   align = "v")
p_all
dev.off()

pdf("Plots/scatterplots/All_plots_drought_violins_only.pdf", height=4, width = 8)
p_violin <- ggarrange(p3, p4, p5, nrow = 1, ncol = 3, common.legend = F,
                      labels = c("A", "B", "C"), vjust = 1)
p_violin
dev.off()

pdf("Plots/scatterplots/All_plots_drought_legendright_violins_new.pdf", height=10, width = 9)
p_all <- ggarrange(p_scat, p_rec, p_violin, ncol = 1, nrow = 3, widths = c(1.25, 1.25, 1), 
                   align = "v")
p_all
dev.off()

##Add arrows

BBCH_wide_last_three <- BBCH_wide[, c(1, 3:4, 6:10)]


pdf("Plots/scatterplots/All_drought_BBCH_30_23_18_cols_arrows.pdf", height=6, width = 8)
p1 <- ggplot(BBCH_wide_last_three, aes(x=BBCH.date30.08, shape =Treat)) + 
  geom_point(aes(y = BBCH.date23.08, colour = "23/08/2023"), size = 2) + 
  geom_point(aes(y = BBCH.date18.08, colour = "30/08/2023"), size = 2) + 
  
  geom_segment(
    aes(
      x = BBCH.date30.08, y = BBCH.date18.08,  # Start point (18/08)
      xend = BBCH.date30.08, yend = BBCH.date23.08,  # End point (23/08)
      #colour = Treat
    ),
    arrow = arrow(length = unit(0.2, "cm")),  # Adds arrowhead
    linewidth = 0.5  # Arrow thickness
  ) +
  
  xlab("BBCH 30/08/2023") + ylab("BBCH") + 
  theme_bw() +
  scale_fill_brewer(palette = "Dark2") +
  # remove space between plot area and x axis
  theme(
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size = 12)) +
  scale_colour_discrete(na.translate = FALSE, name = "Date")  +
  scale_shape_discrete(name = "Treatment")
  #scale_y_continuous(limits = c(20, 90))
p1
dev.off()


######### Box plots by week ---------------
write.xlsx(bbch_all_result, "bbch_All_results_droughtonly.xlsx")

library(ggpubr)
library(dplyr)
library(ggsignif)
library(RColorBrewer)
library(paletteer)
library(ggthemes)
library(rstatix)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcompView)

bbch_all_result$Date <- bbch_all_result$date
bbch_all_result$Date <- gsub('date10.08', 'Week 2', bbch_all_result$Date)
bbch_all_result$Date <- gsub('date30.08', 'Week 5', bbch_all_result$Date)
bbch_all_result$Date <- gsub('date28.07', 'Week 1', bbch_all_result$Date)
bbch_all_result$Date <- gsub('date23.08', 'Week 4', bbch_all_result$Date)
bbch_all_result$Date <- gsub('date18.08', 'Week 3', bbch_all_result$Date)

bbch_all_result$GH <- factor(bbch_all_result$GH)

# mean_FM <- aggregate(. ~ GH + Date,
#                      data = all_FM,
#                      FUN = mean,
#                      na.rm = TRUE)


#Fm'
all_FM <- bbch_all_result[, c(25, 4, 15)]

model <- aov(Fm. ~ GH * Date, data = all_FM)
emm <- emmeans(model, ~ GH | Date)   # sensors within each week
pairs_emm <- pairs(emm, adjust = "holm")
letters <- multcomp::cld(emm, adjust = "holm", Letters = LETTERS, alpha = 0.05)
emm_means <- as.data.frame(emm)
letters_df <- as.data.frame(letters) %>%
  mutate(
    week_num = as.numeric(gsub("[^0-9]", "", as.character(Date))),  # keep only digits
    .group = paste0(.group, "^{", week_num, "}")
  )


Fm <- ggplot(all_FM, aes(x = GH, y = Fm.)) +
  #geom_violin(alpha = 1, width = 3) +
  geom_boxplot(width = 0.6, alpha = 1) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~  Date, ncol = 5, nrow = 1) +
  coord_cartesian(ylim = c(50, 500)) +
  scale_x_discrete(labels = c(
    "Determinate_bush"       = "Drought\nDeterminate Bush", 
    "Indeterminate_bush"     = "Drought\nIndeterminate Bush", 
    "Indeterminate_climbing" = "Drought\nIndeterminate Climbing"
  )) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Stress Condition and Growth Habit", y = "Fm'") +
  geom_text(
    data = letters_df,
    aes(x = GH, y = 500, label = .group),
    inherit.aes = FALSE,
    colour = "black", size = 4,
    parse = TRUE    # tells ggplot to interpret plotmath
  )


pdf("Plots/FM_box_growthhabit_stats_week_ggplot_ANOVA.pdf", height=6, width = 18)
Fm
dev.off()

#Fs

all_FS <- bbch_all_result[, c(25, 4, 11)]

model <- aov(Fs ~ GH * Date, data = all_FS)
emm <- emmeans(model, ~ GH | Date)   # sensors within each week
pairs_emm <- pairs(emm, adjust = "holm")
letters <- multcomp::cld(emm, adjust = "holm", Letters = LETTERS, alpha = 0.05)
emm_means <- as.data.frame(emm)
letters_df <- as.data.frame(letters) %>%
  mutate(
    week_num = as.numeric(gsub("[^0-9]", "", as.character(Date))),  # keep only digits
    .group = paste0(.group, "^{", week_num, "}")
  )

FS <- ggplot(all_FS, aes(x = GH, y = Fs)) +
  #geom_violin(alpha = 1, width = 3) +
  geom_boxplot(width = 0.6, alpha = 1) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~  Date, ncol = 5, nrow = 1) +
  coord_cartesian(ylim = c(50, 260)) +
  scale_x_discrete(labels = c(
    "Determinate_bush"       = "Drought\nDeterminate Bush", 
    "Indeterminate_bush"     = "Drought\nIndeterminate Bush", 
    "Indeterminate_climbing" = "Drought\nIndeterminate Climbing"
  )) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(), 
        strip.text = element_blank()) +
  labs(x = "Stress Condition and Growth Habit", y = "Fs") +
  geom_text(
    data = letters_df,
    aes(x = GH, y = 240, label = .group),
    inherit.aes = FALSE,
    colour = "black", size = 4,
    parse = TRUE    # tells ggplot to interpret plotmath
  )

pdf("Plots/FS_box_growthhabit_stats_week_ggplot_ANOVA.pdf", height=6, width = 18)
FS
dev.off()

#Tleaf
all_TL <- bbch_all_result[, c(25, 4, 10)]

model <- aov(Tleaf ~ GH * Date, data = all_TL)
emm <- emmeans(model, ~ GH | Date)   # sensors within each week
pairs_emm <- pairs(emm, adjust = "holm")
letters <- multcomp::cld(emm, adjust = "holm", Letters = LETTERS, alpha = 0.05)
emm_means <- as.data.frame(emm)
letters_df <- as.data.frame(letters) %>%
  mutate(
    week_num = as.numeric(gsub("[^0-9]", "", as.character(Date))),  # keep only digits
    .group = paste0(.group, "^{", week_num, "}")
  )

TL <- ggplot(all_TL, aes(x = GH, y = Tleaf)) +
  #geom_violin(alpha = 1, width = 3) +
  geom_boxplot(width = 0.6, alpha = 1) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~  Date, ncol = 5, nrow = 1) +
  coord_cartesian(ylim = c(17, 40)) +
  scale_x_discrete(labels = c(
    "Determinate_bush"       = "Drought\nDeterminate Bush", 
    "Indeterminate_bush"     = "Drought\nIndeterminate Bush", 
    "Indeterminate_climbing" = "Drought\nIndeterminate Climbing"
  )) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(), 
        strip.text = element_blank()) +
  labs(x = "Stress Condition and Growth Habit", y = "Tleaf (°C)") +
  geom_text(
    data = letters_df,
    aes(x = GH, y = 39, label = .group),
    inherit.aes = FALSE,
    colour = "black", size = 4,
    parse = TRUE    # tells ggplot to interpret plotmath
  )

pdf("Plots/TL_box_growthhabit_stats_week_ggplot_ANOVA.pdf", height=6, width = 18)
TL
dev.off()

#GSW

all_GS <- bbch_all_result[, c(25, 4, 13)]

model <- aov(gsw ~ GH * Date, data = all_GS)
emm <- emmeans(model, ~ GH | Date)   # sensors within each week
pairs_emm <- pairs(emm, adjust = "holm")
letters <- multcomp::cld(emm, adjust = "holm", Letters = LETTERS, alpha = 0.05)
emm_means <- as.data.frame(emm)
letters_df <- as.data.frame(letters) %>%
  mutate(
    week_num = as.numeric(gsub("[^0-9]", "", as.character(Date))),  # keep only digits
    .group = paste0(.group, "^{", week_num, "}")
  )

GS <- ggplot(all_GS, aes(x = GH, y = gsw)) +
  #geom_violin(alpha = 1, width = 3) +
  geom_boxplot(width = 0.6, alpha = 1) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~  Date, ncol = 5, nrow = 1) +
  coord_cartesian(ylim = c(0, 1.2)) +
  scale_x_discrete(labels = c(
    "Determinate_bush"       = "Drought\nDeterminate Bush", 
    "Indeterminate_bush"     = "Drought\nIndeterminate Bush", 
    "Indeterminate_climbing" = "Drought\nIndeterminate Climbing"
  )) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Stress Condition and Growth Habit", y = expression(Gsw ~ (mol ~ m^{-2} ~ s^{-1}))) +
  geom_text(
    data = letters_df,
    aes(x = GH, y = 1.1, label = .group),
    inherit.aes = FALSE,
    colour = "black", size = 4,
    parse = TRUE    # tells ggplot to interpret plotmath
  )

pdf("Plots/GS_box_growthhabit_stats_week_ggplot_ANOVA.pdf", height=6, width = 18)
GS
dev.off()

#E apparent 

all_E <- bbch_all_result[, c(25, 4, 14)]

model <- aov(E_apparent ~ GH * Date, data = all_E)
emm <- emmeans(model, ~ GH | Date)   # sensors within each week
pairs_emm <- pairs(emm, adjust = "holm")
letters <- multcomp::cld(emm, adjust = "holm", Letters = LETTERS, alpha = 0.05)
emm_means <- as.data.frame(emm)
letters_df <- as.data.frame(letters) %>%
  mutate(
    week_num = as.numeric(gsub("[^0-9]", "", as.character(Date))),  # keep only digits
    .group = paste0(.group, "^{", week_num, "}")
  )

EA <- ggplot(all_E, aes(x = GH, y = E_apparent)) +
  #geom_violin(alpha = 1, width = 3) +
  geom_boxplot(width = 0.6, alpha = 1) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~  Date, ncol = 5, nrow = 1) +
  coord_cartesian(ylim = c(-0.3, 7)) +
  scale_x_discrete(labels = c(
    "Determinate_bush"       = "Drought\nDeterminate Bush", 
    "Indeterminate_bush"     = "Drought\nIndeterminate Bush", 
    "Indeterminate_climbing" = "Drought\nIndeterminate Climbing"
  )) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(), 
        strip.text = element_blank()) +
  labs(x = "Stress Condition and Growth Habit", y = expression(italic(E) * " apparent" ~ (mmol ~ m^{-2} ~ s^{-1}))) +
  geom_text(
    data = letters_df,
    aes(x = GH, y = 7, label = .group),
    inherit.aes = FALSE,
    colour = "black", size = 4,
    parse = TRUE    # tells ggplot to interpret plotmath
  )

pdf("Plots/E_box_growthhabit_stats_week_ggplot_ANOVA.pdf", height=6, width = 18)
EA
dev.off()

#PhiPS2 

all_PI <- bbch_all_result[, c(25, 4, 16)]

model <- aov(PhiPS2 ~ GH * Date, data = all_PI)
emm <- emmeans(model, ~ GH | Date)   # sensors within each week
pairs_emm <- pairs(emm, adjust = "holm")
letters <- multcomp::cld(emm, adjust = "holm", Letters = LETTERS, alpha = 0.05)
emm_means <- as.data.frame(emm)
letters_df <- as.data.frame(letters) %>%
  mutate(
    week_num = as.numeric(gsub("[^0-9]", "", as.character(Date))),  # keep only digits
    .group = paste0(.group, "^{", week_num, "}")
  )

PI <- ggplot(all_PI, aes(x = GH, y = PhiPS2)) +
  #geom_violin(alpha = 1, width = 3) +
  geom_boxplot(width = 0.6, alpha = 1) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~  Date, ncol = 5, nrow = 1) +
  coord_cartesian(ylim = c(0, 0.8)) +
  scale_x_discrete(labels = c(
    "Determinate_bush"       = "Drought\nDeterminate Bush", 
    "Indeterminate_bush"     = "Drought\nIndeterminate Bush", 
    "Indeterminate_climbing" = "Drought\nIndeterminate Climbing"
  )) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(), 
        strip.text = element_blank()) +
  labs(x = "Stress Condition and Growth Habit", y = expression(PhiPS2 ~ (mu * mol ~ m^{-2} ~ s^{-1}))) +
  geom_text(
    data = letters_df,
    aes(x = GH, y = 0.8, label = .group),
    inherit.aes = FALSE,
    colour = "black", size = 4,
    parse = TRUE    # tells ggplot to interpret plotmath
  )

pdf("Plots/PI_box_growthhabit_stats_week_ggplot_ANOVA.pdf", height=6, width = 18)
PI
dev.off()


pdf("Plots/allgrowthhabit_bar_ggplot_weekA.pdf", height=12, width = 16) #originally 20 and 15 portrait, 18 and 20 for landscape
p2 <- ggarrange(Fm, FS, TL, 
                ncol = 1, 
                nrow = 3, 
                legend = "right",
                labels = "AUTO"
)
p2
dev.off()

pdf("Plots/allgrowthhabit_bar_ggplot_weekB.pdf", height=12, width = 16) #originally 20 and 15 portrait, 18 and 20 for landscape
p2 <- ggarrange(GS, EA, PI, 
                ncol = 1, 
                nrow = 3, 
                legend = "right",
                labels = c("D", "E", "F")
)
p2
dev.off()


#ETR
all_etr <- bbch_all_result[, c(25, 4, 8)]

model <- aov(ETR ~ GH * Date, data = all_etr)
emm <- emmeans(model, ~ GH | Date)   # sensors within each week
pairs_emm <- pairs(emm, adjust = "holm")
letters <- multcomp::cld(emm, adjust = "holm", Letters = LETTERS, alpha = 0.05)
emm_means <- as.data.frame(emm)
letters_df <- as.data.frame(letters) %>%
  mutate(
    week_num = as.numeric(gsub("[^0-9]", "", as.character(Date))),  # keep only digits
    .group = paste0(.group, "^{", week_num, "}")
  )

