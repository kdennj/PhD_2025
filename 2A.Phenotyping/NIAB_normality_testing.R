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

setwd("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Bean project/23.03-Field-trial-NIAB/Licor/")

###Normal histograms -------------------
date30.08 <- read.xlsx("30.08.23.merged.xlsx")
date23.08 <- read.xlsx("23.08.23.merged.xlsx")
date18.08 <- read.xlsx("18.08.23.merged.xlsx")
date10.08 <- read.xlsx("10.08.23.merged.xlsx")
date28.07 <- read.xlsx("28.07.23.merged.xlsx")


df_list <- list(date30.08, date23.08, date18.08, date10.08, date28.07)
names(df_list) <- c("date30.08", "date23.08", "date18.08", "date10.08", "date28.07")
column_selection <- c(111, 37, 23, 41, 32, 38, 16, 19, 33, 34)

# Loop through each dataframe in the list
for (i in seq_along(df_list)) {
  # Add the new 'Accession_T' column
  df_list[[i]]$Accession_T <- paste(df_list[[i]]$Accession, df_list[[i]]$Treat, sep = "_")
  
  # Select and reorder the specified columns (if they exist)
  valid_columns <- column_selection[column_selection %in% seq_len(ncol(df_list[[i]]))]
  df_list[[i]] <- df_list[[i]][, c(valid_columns, ncol(df_list[[i]])), drop = FALSE]
  
  #Remove last column 
  df_list[[i]] <- df_list[[i]][, -ncol(df_list[[i]]), drop = FALSE]
  
  # Add the 'date' column using the dataframe name
  df_list[[i]]$date <- names(df_list)[i]
  
  
}

#which(colnames(df_list[[1]]) == "30/08/2023_BBCH")

all_df <- do.call(rbind, df_list)
rownames(all_df) <- NULL
write.xlsx(all_df, "Plots/corrplots/all_dates_licor.xlsx")

all_df[, -c(1, 11)] <- lapply(all_df[, -c(1, 11)], as.numeric)
all_df <- all_df[, c(-11)]

result <- aggregate(. ~ Accession_T, 
                    data = all_df, 
                    FUN = mean, 
                    na.rm = TRUE)

result <- result[!grepl("_C", result$Accession_T), ]


write.xlsx(result, "Plots/corrplots/all_dates_licor_mean_df_drought.xlsx")

#Plotting
output_dir <- "Plots/corrplots/histograms"


# pdf("Plots/corrplots/histograms/ETR_hist.pdf", height=6, width = 14)
# p1 <- ggplot(ETR_mean, aes(x=mean)) + geom_histogram(bins = 10) + 
#     xlab("gsw") + ylab("Count") + theme_bw() + 
#   # remove space between plot area and x axis
#   theme(
#     axis.title.y = element_text(size=12),
#     axis.title.x = element_text(size = 12))
# p1
# dev.off()


for (col_name in colnames(all_df)[-1]) {
  # Create the plot
  p <- ggplot(all_df, aes_string(x = col_name)) +
    geom_histogram(bins = 10, fill = "blue", color = "black") +
    xlab(col_name) +
    ylab("Count") +
    theme_bw() +
    theme(
      axis.title.y = element_text(size = 12),
      axis.title.x = element_text(size = 12)
    )
  
  # Save the plot to a file named after the column
  ggsave(
    filename = paste0(output_dir, "/", col_name, "_histogram.pdf"),
    plot = p,
    height = 6,
    width = 14
  )
}

#### Shapiro wilk test for normality ----
library("ggpubr")

#qq plots
output_dir <- "Plots/corrplots/histograms"

for (col_name in colnames(all_df)[-1]) {
  # Create the plot
  p <- ggqqplot(all_df[[col_name]]) + theme_bw()
  # Save the plot to a file named after the column
  ggsave(
    filename = paste0(output_dir, "/", col_name, "_qqplot.pdf"),
    plot = p,
    height = 6,
    width = 6
  )
}

#Shapiro wilk test

for (col_name in colnames(all_df)[-1]) {
  # Perform Shapiro-Wilk test
  p <- shapiro_test((all_df[[col_name]])) 
  
  # Check if the p-value is greater than 0.05
  if (p$p.value > 0.05) {
    print(paste("Column:", col_name, "- p-value:", p$p.value))
  }
}

# Loop through each column (except the first)
for (col_name in colnames(all_df)[-1]) {
  # Fit a simple linear model with the column as the response and an intercept-only model
  model <- lm(all_df[[col_name]] ~ 1)
  
  # Extract residuals
  residuals <- resid(model)
  
  # Perform Shapiro-Wilk test on residuals
  p <- shapiro.test(residuals)
  
  # Check if the p-value is greater than 0.05
  if (p$p.value > 0.05) {
    print(paste("Column:", col_name, "- p-value:", p$p.value))
  }
}

###Add weights for normality -----------------

harvest <- read_excel("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Bean project/23.03-Field-trial-NIAB/KDJ Bean panel randomisation.xlsx", sheet = "harvest_edit")
colnames(harvest)[7] <- c("Pot.Position.GH")
harvest <- harvest[, c(5, 10, 12:14)]

harvest$Accession_T <- paste(harvest$Accession, harvest$Treat, sep = "_")
harvest <- harvest[, c(6, 3, 4, 5)]

harvest[, -c(1)] <- lapply(harvest[, -c(1)], as.numeric)

result_weights <- aggregate(. ~ Accession_T, 
                            data = harvest, 
                            FUN = mean, 
                            na.rm = TRUE)

write.xlsx(result_weights, "Plots/corrplots/result_weights_mean_all.xlsx")

result_weights <- result_weights[!grepl("_C", result_weights$Accession_T), ]
colnames(result_weights) <- c("Accession_T", "No_pods", "Pod_weight", "All_foliar_weight")

output_dir <- "Plots/corrplots/histograms"

for (col_name in colnames(result_weights)[-1]) {
  # Create the plot
  p <- ggplot(result_weights, aes_string(x = col_name)) +
    geom_histogram(bins = 10, fill = "blue", color = "black") +
    xlab(col_name) +
    ylab("Count") +
    theme_bw() +
    theme(
      axis.title.y = element_text(size = 12),
      axis.title.x = element_text(size = 12)
    )
  
  # Save the plot to a file named after the column
  ggsave(
    filename = paste0(output_dir, "/", col_name, "_histogram.pdf"),
    plot = p,
    height = 6,
    width = 14
  )
}

for (col_name in colnames(result_weights[-1])) {
  # Create the plot
  p <- ggqqplot(result_weights[[col_name]]) + theme_bw()
  # Save the plot to a file named after the column
  ggsave(
    filename = paste0(output_dir, "/", col_name, "_qqplot.pdf"),
    plot = p,
    height = 6,
    width = 6
  )
}

# Loop through each column (except the first)
for (col_name in colnames(result_weights)[-1]) {
  # Fit a simple linear model with the column as the response and an intercept-only model
  model <- lm(result_weights[[col_name]] ~ 1)
  
  # Extract residuals
  residuals <- resid(model)
  
  # Perform Shapiro-Wilk test on residuals
  p <- shapiro.test(residuals)
  
  # Check if the p-value is greater than 0.05
  if (p$p.value > 0.05) {
    print(paste("Column:", col_name, "- p-value:", p$p.value))
  }
}

