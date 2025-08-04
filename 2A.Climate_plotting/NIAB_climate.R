library("readxl")
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)
library(plyr)
library(tidyverse)
library(data.table)
library(purrr)
library(ggstatsplot)
library(ggplot2)
library(ggrepel)
library(ggsignif)

#Thermal time 
#Base temperature common beans 12 degrees
#Real temperature for the days - average of the days with BBCH data 

setwd("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Bean project/23.03-Field-trial-NIAB/Visual_crossing/")

temp_data <- read.xlsx("cambridge 2023-06-19 to 2023-09-09.xlsx")

precip <- temp_data[, c(2, 11)]
all <- temp_data[, c(2, 5, 11, 27, 28)] #10 is humidity 
all$sunset <- gsub("T", " ", all$sunset)
all$sunrise <- gsub("T", " ", all$sunrise)

all$sunrise <- as.POSIXct(all$sunrise, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
all$sunset <- as.POSIXct(all$sunset, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Calculate the day length
all$daylength_sec <- as.numeric(difftime(all$sunset, all$sunrise, units = "secs"))
all$daylength_hours <- all$daylength_sec / 3600

all <- all[, c(1:3, 7)]

long <- melt(setDT(all), id.vars = c("datetime"), variable.name = "climate")

pdf("../Plots/precip.pdf", height=7, width = 10)
p <- ggplot(precip, aes(x=datetime, y=precip, group = 1)) +
  geom_line() +
  geom_point()+
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
dev.off()

xlabs <- long$datetime[seq(1, nrow(long), by = 2)]

x <- long$datetime[39]

str(long)
long$datetime <- as.Date(long$datetime)

# Determine y-axis limits without gaps
ylim <- range(temps)
ylim <- c(floor(min(ylim)), ceiling(max(ylim)))

pdf("../Plots/precip_temp_hum_fewer_lines.pdf", height=5, width = 12)
p <- ggplot(long, aes(x=datetime, y=value, group = climate, colour = climate)) +
  geom_line() +
  geom_point()+
  theme_bw() + scale_color_brewer(palette = "Dark2", labels = c("Temperature", "Humidity", "Precipitation")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_vline(xintercept =as.Date("2023-07-27"), linetype="dotted", linewidth = 0.5) +
  geom_text(aes(x=as.Date("2023-07-27"), label="\nIrrigation stopped", y=50), colour="black", angle=90) +
  geom_vline(xintercept =as.Date("2023-08-25"), linetype="dotted", linewidth = 0.5) +
  geom_text(aes(x=as.Date("2023-08-25"), label="\nRecovery", y=50), colour="black", angle=90) +
  labs(x = "Date", y = "Value") 
p
dev.off()


##Plotting precip, temp, daylength 
str(long)

colnames(long)[2] <- c("Climate")

pdf("../Plots/precip_temp_daylength_fewer_lines2.pdf", height=5, width = 8)
p <- ggplot(long, aes(x=datetime, y=value, group = Climate, colour = Climate)) +
  geom_line() +
  geom_point()+
  theme_bw() + scale_color_brewer(palette = "Dark2", labels = c("Temperature (°C)", "Precipitation (mm)", "Daylength (hrs)")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_vline(xintercept =as.Date("2023-07-27"), linetype="dotted", linewidth = 1) +
  #geom_text(aes(x=as.Date("2023-07-27"), label="\nIrrigation stopped", y=25, family = "sans", fontface = "plain"), colour="black", angle=0, hjust = 1, size = 4) +
  annotate("text",   x = as.Date("2023-07-27"), y = 25, label = "\nIrrigation stopped", colour = "black",
           angle = 0, hjust = 1, size = 4, family = "sans", fontface = "plain") +
  geom_vline(xintercept =as.Date("2023-08-25"), linetype="dotted", linewidth = 1) +
  #geom_text(aes(x=as.Date("2023-08-25"), label="\nRecovery", y=25), colour="black", angle=0, hjust = 1, size = 4) +
  annotate("text",   x = as.Date("2023-08-25"), y = 25, label="\nRecovery", colour = "black",
           angle = 0, hjust = 1, size = 4, family = "sans", fontface = "plain") +
  labs(x = "Date", y = "Value") 
p
dev.off()

##Add vertical lines for removing water / recovery rain? 

#27/07/2023 water removed
# 25/08/2023 recovery

#Calculate mean values
mean(all$temp)
mean(all$precip)
mean(all$daylength_hours)

irrigation_stopped <- all[all$datetime >= "2023-07-27", ]
recovery <- all[all$datetime >= "2023-08-25", ]

# Calculate the mean temperature
mean(recovery$precip)

#### Sensors -------------------------------------

#Import the data frames
Sensor_1 <- read.xlsx("../Sensors/sensor1/sensor1.xlsx", sheet = "Processed Data Config 3")
Sensor_2 <- read.xlsx("../Sensors/sensor2/sensor2.xlsx", sheet = "Processed Data Config 4")
Sensor_3 <- read.xlsx("../Sensors/sensor3/sensor3B.xlsx", sheet = "Processed Data Config 3")
Sensor_4 <- read.xlsx("../Sensors/sensor4/sensor4B.xlsx", sheet = "Processed Data Config 2")
Sensor_5 <- read.xlsx("../Sensors/sensor5/sensor5.xlsx", sheet = "Processed Data Config 3") #Removed one of the Matric Potential due to large values
Sensor_6 <- read.xlsx("../Sensors/sensor6/sensor6.xlsx", sheet = "Processed Data Config 3")

#Make function to convert timestamp column to date and time
convert_sensor_dates <- function(df) {
  for (col in names(df)) {
    if (grepl("Timestamp", col, ignore.case = TRUE) && is.numeric(df[[col]])) {
      df[[col]] <- as.POSIXct(df[[col]] * 86400, origin = "1899-12-30")
    }
  }
  return(df)
}

# Loop through all data frames with 'sensor' in the name
for (df_name in ls()) {
  if (grepl("sensor", df_name, ignore.case = TRUE)) {
    # Get the data frame using get()
    df <- get(df_name)
    # Apply the conversion function
    df <- convert_sensor_dates(df)
    # Assign it back to the same name
    assign(df_name, df)
  }
}

#Select growth habit and stress
Control_determinate <- Sensor_1[, c(-14, -15)]
Sens2_drought_determinate <- Sensor_2[, c(-14, -15)]
Sens3_drought_determinate <- Sensor_3[, c(1:6)]
Drought_determinate <- merge(Sens2_drought_determinate, Sens3_drought_determinate, by = "Timestamp", all = TRUE)
Sens3_drought_indeterminatebush <- Sensor_3[, c(1, 7:11)]
Sens4_drought_indeterminatebush <- Sensor_4[, c(1:6)]
Drought_indeterminatebush <- merge(Sens3_drought_indeterminatebush, Sens4_drought_indeterminatebush, by = "Timestamp", all = TRUE)
Sens4_drought_indeterminateclimb <- Sensor_4[ , c(1, 7:13)]
Sens5_drought_indeterminateclimb <- Sensor_5[, c(1:10)]
Sens6_drought_indeterminateclimb <- Sensor_6[, c(1:11)]
D_IC <- list(Sens4_drought_indeterminateclimb, Sens5_drought_indeterminateclimb, Sens6_drought_indeterminateclimb)
Drought_indeterminateclimb <- Reduce(function(x, y) merge(x, y, by = "Timestamp", all = TRUE), D_IC)

write.xlsx(Control_determinate, "../Plots/Sensors/Control_determinate_sensors.xlsx")
write.xlsx(Drought_determinate, "../Plots/Sensors/Drought_determinate_sensors.xlsx")
write.xlsx(Drought_indeterminatebush, "../Plots/Sensors/Drought_indeterminatebush_sensors.xlsx")
write.xlsx(Drought_indeterminateclimb, "../Plots/Sensors/Drought_indeterminateclimbing_sensors.xlsx")

## Calculate means for the different measurements
#Function to create average columns
variables <- c("Matric.Potential", "Soil.Temperature", "Bulk.EC", "Water.Content")

calculate_averages <- function(data, variables, df_name) {
  for (var in variables) {
    avg_col_name <- paste0(var, "_avg_", df_name)  # Append dataframe name
    data[[avg_col_name]] <- data %>%
      select(matches(var)) %>%
      rowMeans(na.rm = TRUE)
  }
  return(data)
}

# # Loop through Sensor_1 to Sensor_6
# for (i in 1:6) {
#   sensor_name <- paste0("Sensor_", i)  # Create the dataframe name
#   df <- get(sensor_name)                # Get the dataframe
#   
#   # Calculate averages and append the dataframe name
#   df <- calculate_averages(df, variables, sensor_name)
#   
#   # Assign the modified dataframe back to the original variable name
#   assign(sensor_name, df)
# }

#Loop throught for dataframes by growth habit and stress
dfs <- list(
  Control_determinate = Control_determinate, 
  Drought_determinate = Drought_determinate, 
  Drought_indeterminatebush = Drought_indeterminatebush, 
  Drought_indeterminateclimb = Drought_indeterminateclimb
)
#List to store updated data frames
updated_dfs <- list()

# Iterate over each data frame in the list
for (sensor_name in names(dfs)) {
  df <- dfs[[sensor_name]]
  df <- calculate_averages(df, variables, sensor_name)
  updated_dfs[[sensor_name]] <- df  # Store in a new list
}

# Access the updated data frames
Control_determinate <- updated_dfs$Control_determinate  
Drought_determinate <- updated_dfs$Drought_determinate
Drought_indeterminatebush <- updated_dfs$Drought_indeterminatebush
Drought_indeterminateclimb <- updated_dfs$Drought_indeterminateclimb

##Combine all into 1; Matric Potential --------
# MP_1 <- Sensor_1[, c(1, 16)]
# MP_2 <- Sensor_2[, c(1, 16)]
# MP_3 <- Sensor_3[, c(1, 14)]
# MP_4 <- Sensor_4[, c(1, 16)]
# MP_5 <- Sensor_5[, c(1, 15)]
# MP_6 <- Sensor_6[, c(1, 16)]

MP_Control_determinate <- Control_determinate[, c(1, 14)]
MP_Drought_determinate <- Drought_determinate[, c(1, 19)]
MP_Drought_indeterminatebush <- Drought_indeterminatebush[, c(1, 12)]
MP_Drought_indeterminateclimb <- Drought_indeterminateclimb[, c(1, 28)]

dfs <- list(MP_Control_determinate, MP_Drought_determinate, MP_Drought_indeterminatebush, MP_Drought_indeterminateclimb) 
MP_all <- Reduce(function(x, y) merge(x, y, by = "Timestamp", all = TRUE), dfs)

#colnames(MP_all) <- c("Timestamp", "Sensor_1_C", "Sensor_2_D", "Sensor_3_D", "Sensor_4_D", "Sensor_5_D", "Sensor_6_D")
colnames(MP_all) <- c("Timestamp", "Control_determinatebush", "Drought_determinatebush", "Drought_indeterminatebush", "Drought_indeterminateclimbing")


##Combine all into 1; Soil temperature --------
# ST_1 <- Sensor_1[, c(1, 17)]
# ST_2 <- Sensor_2[, c(1, 17)]
# ST_3 <- Sensor_3[, c(1, 15)]
# ST_4 <- Sensor_4[, c(1, 17)]
# ST_5 <- Sensor_5[, c(1, 16)]
# ST_6 <- Sensor_6[, c(1, 17)]

ST_Control_determinate <- Control_determinate[, c(1, 15)]
ST_Drought_determinate <- Drought_determinate[, c(1, 20)]
ST_Drought_indeterminatebush <- Drought_indeterminatebush[, c(1, 13)]
ST_Drought_indeterminateclimb <- Drought_indeterminateclimb[, c(1, 29)]

dfs <- list(ST_Control_determinate, ST_Drought_determinate, ST_Drought_indeterminatebush, ST_Drought_indeterminateclimb) 
ST_all <- Reduce(function(x, y) merge(x, y, by = "Timestamp", all = TRUE), dfs)
colnames(ST_all) <- c("Timestamp", "Control_determinatebush", "Drought_determinatebush", "Drought_indeterminatebush", "Drought_indeterminateclimbing")


# dfs <- list(ST_1, ST_2, ST_3, ST_4, ST_5, ST_6) 
# st_all <- Reduce(function(x, y) merge(x, y, by = "Timestamp", all = TRUE), dfs)

##Combine all into 1; Bulk EC --------
# EC_1 <- Sensor_1[, c(1, 18)]
# EC_2 <- Sensor_2[, c(1, 18)]
# EC_3 <- Sensor_3[, c(1, 16)]
# EC_4 <- Sensor_4[, c(1, 18)]
# EC_5 <- Sensor_5[, c(1, 17)]
# EC_6 <- Sensor_6[, c(1, 18)]

EC_Control_determinate <- Control_determinate[, c(1, 16)]
EC_Drought_determinate <- Drought_determinate[, c(1, 21)]
EC_Drought_indeterminatebush <- Drought_indeterminatebush[, c(1, 14)]
EC_Drought_indeterminateclimb <- Drought_indeterminateclimb[, c(1, 30)]

dfs <- list(EC_Control_determinate, EC_Drought_determinate, EC_Drought_indeterminatebush, EC_Drought_indeterminateclimb) 
EC_all <- Reduce(function(x, y) merge(x, y, by = "Timestamp", all = TRUE), dfs)
colnames(EC_all) <- c("Timestamp", "Control_determinatebush", "Drought_determinatebush", "Drought_indeterminatebush", "Drought_indeterminateclimbing")


##Combine all into 1; Water content --------
# WC_1 <- Sensor_1[, c(1, 19)]
# WC_2 <- Sensor_2[, c(1, 19)]
# WC_3 <- Sensor_3[, c(1, 17)]
# WC_4 <- Sensor_4[, c(1, 19)]
# WC_5 <- Sensor_5[, c(1, 18)]
# WC_6 <- Sensor_6[, c(1, 19)]

WC_Control_determinate <- Control_determinate[, c(1, 17)]
WC_Drought_determinate <- Drought_determinate[, c(1, 22)]
WC_Drought_indeterminatebush <- Drought_indeterminatebush[, c(1, 15)]
WC_Drought_indeterminateclimb <- Drought_indeterminateclimb[, c(1, 31)]

dfs <- list(WC_Control_determinate, WC_Drought_determinate, WC_Drought_indeterminatebush, WC_Drought_indeterminateclimb) 
WC_all <- Reduce(function(x, y) merge(x, y, by = "Timestamp", all = TRUE), dfs)
colnames(WC_all) <- c("Timestamp", "Control_determinatebush", "Drought_determinatebush", "Drought_indeterminatebush", "Drought_indeterminateclimbing")


#Change into long format
long_MP <- melt(setDT(MP_all), id.vars = c("Timestamp"), variable.name = "Sensor")
long_ST <- melt(setDT(ST_all), id.vars = c("Timestamp"), variable.name = "Sensor")
long_EC <- melt(setDT(EC_all), id.vars = c("Timestamp"), variable.name = "Sensor")
long_WC <- melt(setDT(WC_all), id.vars = c("Timestamp"), variable.name = "Sensor")

write.xlsx(long, "../Plots/Sensors/WC_all_growthhabit.xlsx")

######### Can start from here for sensors and plotting ----------
#Read in excels 
long <- read.xlsx("../Plots/Sensors/MP_all.xlsx")


#### Plotting ----------------------------------
start_date <- as.POSIXct("2023-07-06 00:00:00")
end_date <- as.POSIXct("2023-09-09 00:00:00")

long_filtered_MP <- long_MP[as.Date(long_MP$Timestamp) != as.Date("2023-09-08"), ]
long_filtered_MP <- long_filtered_MP[as.Date(long_filtered_MP$Timestamp) != as.Date("2023-09-07"), ]
long_filtered_MP <- long_filtered_MP[as.Date(long_filtered_MP$Timestamp) != as.Date("2023-09-06"), ]

#Plotting Matric Potential
pdf("../Plots/Sensors/Matric_potential_all_growthhabit.pdf", height=5, width = 12)
p_MP <- ggplot(long_filtered_MP, aes(x=Timestamp, y=value, group = Sensor, colour = Sensor)) +
  geom_line() +
  geom_point(size = 1) + ylab("Soil Matric Potential (kPa)") + xlab("Date") +
  theme_bw() + scale_color_brewer(palette = "Dark2", labels = c("Control determinate bush", "Drought determinate bush", "Drought indeterminate bush", "Drought indeterminate climbing")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_vline(xintercept = as.POSIXct("2023-07-27 00:00:00"), linetype="dotted", linewidth = 1) +
  geom_text(aes(x=as.POSIXct("2023-07-27 00:00:00"), label="\nIrrigation stopped", y=-8000), colour="black", angle=0, hjust = 1) +
  geom_vline(xintercept = as.POSIXct("2023-08-25 00:00:00"), linetype="dotted", linewidth = 1) +
  geom_text(aes(x= as.POSIXct("2023-08-25 00:00:00"), label="\nRecovery", y=-8000), colour="black", angle=0, hjust = 1) +
  scale_x_datetime(breaks = "3 days")
  #+ scale_x_continuous(expand = expansion(mult = c(0.1, 0.1)))
p_MP
dev.off()


#Violin plot Matric Potential
long_filtered_MP_more <- long_filtered_MP %>% filter(Timestamp >= "2023-07-27")

pdf("../Plots/Sensors/Matric_Potential_all_growthhabit_stats_irrigstopped.pdf", height=8, width = 8)
ps_MP <- ggstatsplot::ggbetweenstats(data = long_filtered_MP_more,
                                 x = Sensor,
                                 y = value,
                                 type = "parametric",
                                 pairwise.display = "significant",
                                 p.adjust.method = "holm",
                                 results.subtitle = FALSE,
                                 centrality.point.args = list(size = 4, color = "darkred"),
                                 violin.args = list(width = 1, alpha = 0.2, na.rm = TRUE),
                                 ggsignif.args = list(textsize = 3, tip_length = 0.01, na.rm = TRUE),
                                 # mean.ci = TRUE,
                                 pairwise.comparisons = TRUE,
                                 bf.message = FALSE,
                                 xlab = "Sensors", 
                                 ylab = "Soil Matric Potential (kPa)",
                                 digits = "signif",
                                 ggtheme = ggplot2::theme_bw()
                                 ) +
  scale_x_discrete(labels = c(
    "Control_determinatebush" = "Control \nDeterminate Bush", 
    "Drought_determinatebush" = "Drought \nDeterminate Bush", 
    "Drought_indeterminatebush" = "Drought \nIndeterminate Bush", 
    "Drought_indeterminateclimbing" = "Drought \nIndeterminate Climbing"
   )) #+
  # theme(
  #   axis.text.x = element_text(angle = 45, hjust = 1))
ps_MP
dev.off()

##Bar plot
#https://stackoverflow.com/questions/77602454/how-to-use-ggbetweenstats-with-grouped-data-at-different-timepoints-in-r
#https://www.crumplab.com/psyc7709_2019/book/docs/how-to-annotate-a-graph-using-gg-signif.html


#comparisons <- combn(levels(factor(long$Sensor)), 2, simplify = FALSE)

anova_result <- aov(value ~ Sensor, data = long)
summary(anova_result)
tukey_result <- TukeyHSD(anova_result)
tukey_df <- as.data.frame(tukey_result$Sensor)
tukey_df <- tukey_df[tukey_df$`p adj`<0.05,]

annotations <- scales::pvalue(tukey_df$`p adj`, add_p = TRUE)
y_position <- seq(from = 100, by = 400, length.out = 15)

pdf("../Plots/Sensors/Matric_Potential_all_barplot.pdf", height=12, width = 10)
p <- ggplot(long, aes(x = Sensor, y = value)) +
  geom_point(aes(color = Sensor, fill = after_scale(alpha(colour, 0.5))), 
             position = position_jitterdodge(dodge.width = 0.9, 0.1),
             size = 2, shape = 21) + ylab("Soil Matric Potential (kPa)") +
  geom_boxplot(fill = NA, color = "black", width = 0.2, linewidth = 0.4,
               position = position_dodge(0.9)) +
  geom_point(stat = "summary", size = 3, color = "#8a0f00",
             position = position_dodge(0.9), fun = mean)  +
  geom_label_repel(stat = "summary", fun = mean, size = 3,
                   aes(label = paste0("hat(mu)*scriptstyle(mean)==", 
                                      round(after_stat(y), 2))),
                   parse = TRUE, position = position_dodge(0.9)) +
  scale_y_continuous(sec.axis = sec_axis(~., 
                                         bquote(Pairwise~Test~paste(":")~bold(ANOVA)))) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 10) +
  theme(axis.title = element_text(face = 2),
        legend.position = "bottom",
        axis.text.y.right = element_blank())

  
    # geom_signif(y_position = seq(from = 100, by = 400, length.out = 13), 
  #             annotations = cut(
  #             tukey_df$`p adj`,
  #             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  #             labels = c("***", "**", "*", "")),
  #             tip_length = 0.01, 
  #             xmin = 1:13 - 0.15, 
  #             xmax = 1:13 + 0.15)

# geom_signif(y_position = seq(from = 100, by = 400, length.out = 15), 
#             annotations = scales::pvalue(tukey_df$`p adj`, add_p = TRUE), tip_length = 0.01, 
#             xmin = 1:15 - 0.22, xmax = 1:15 + 0.22)

# geom_signif(y_position = seq(from = 100, by = 400, length.out = 13), 
#             annotations = cut(
#               tukey_df$`p adj`,
#               breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
#               labels = c("***", "**", "*", "")),
#             tip_length = 0.01, 
#             xmin = 1:13 - 0.15, xmax = 1:13 + 0.15)

for (i in seq_along(comparisons)) {
  p <- p + geom_signif(
    comparisons = list(comparisons[[i]]),
    map_signif_level = TRUE,
    y_position = max(long$value, na.rm = TRUE) + i * 100,  # Stagger y-positions
    tip_length = 0.01
  )
}

  
p
dev.off()


###Plotting Soil Temperature --------------------------------------------------
long_filtered_ST <- long_ST[as.Date(long_ST$Timestamp) != as.Date("2023-09-08"), ]
long_filtered_ST <- long_filtered_ST[as.Date(long_filtered_ST$Timestamp) != as.Date("2023-09-07"), ]
long_filtered_ST <- long_filtered_ST[as.Date(long_filtered_ST$Timestamp) != as.Date("2023-09-06"), ]


pdf("../Plots/Sensors/Soil_temp_all_growthhabit.pdf", height=5, width = 12)
p_ST <- ggplot(long_filtered_ST, aes(x=Timestamp, y=value, group = Sensor, colour = Sensor)) +
  geom_line() +
  geom_point(size = 1) + ylab("Soil Temperature (°C)") + xlab("Date") +
  theme_bw() + scale_color_brewer(palette = "Dark2", labels = c("Control determinate bush", "Drought determinate bush", "Drought indeterminate bush", "Drought indeterminate climbing")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_vline(xintercept = as.POSIXct("2023-07-27 00:00:00"), linetype="dotted", linewidth = 1) +
  geom_text(aes(x=as.POSIXct("2023-07-27 00:00:00"), label="\nIrrigation stopped", y=40), colour="black", angle=0, hjust = 0) +
  geom_vline(xintercept = as.POSIXct("2023-08-25 00:00:00"), linetype="dotted", linewidth = 1) +
  geom_text(aes(x= as.POSIXct("2023-08-25 00:00:00"), label="\nRecovery", y=40), colour="black", angle=0, hjust = 0) +
  scale_x_datetime(breaks = "3 days")
p_ST
dev.off()

#Violin plot Soil temperature
long_filtered_ST_more <- long_filtered_ST %>% filter(Timestamp >= "2023-07-27")


pdf("../Plots/Sensors/Soil_temp_all_growthhabit_stats_irrigstopped.pdf", height=8, width = 8)
ps_ST <- ggstatsplot::ggbetweenstats(data = long_filtered_ST_more,
                                 x = Sensor,
                                 y = value,
                                 type = "parametric",
                                 pairwise.display = "significant",
                                 pairwise.annotation = "asterisk",
                                 p.adjust.method = "holm",
                                 results.subtitle = FALSE,
                                 centrality.point.args = list(size = 4, color = "darkred"),
                                 violin.args = list(width = 0.5, alpha = 0.2, na.rm = TRUE),
                                 ggsignif.args = list(textsize = 3, tip_length = 0.01, na.rm = TRUE),
                                 # mean.ci = TRUE,
                                 pairwise.comparisons = TRUE,
                                 bf.message = FALSE,
                                 xlab = "Sensors", 
                                 ylab = "Soil Temperature (°C)",
                                 digits = "signif",
                                 ggtheme = ggplot2::theme_bw(), 
                                 boxplot.args = list(width = 0) 
  ) +
  scale_x_discrete(labels = c(
    "Control_determinatebush" = "Control \nDeterminate Bush", 
    "Drought_determinatebush" = "Drought \nDeterminate Bush", 
    "Drought_indeterminatebush" = "Drought \nIndeterminate Bush", 
    "Drought_indeterminateclimbing" = "Drought \nIndeterminate Climbing"))
ps_ST
dev.off()

###Plotting EC ------------------------------------
long_filtered_EC <- long_EC[as.Date(long_EC$Timestamp) != as.Date("2023-09-08"), ]
long_filtered_EC <- long_filtered_EC[as.Date(long_filtered_EC$Timestamp) != as.Date("2023-09-07"), ]
long_filtered_EC <- long_filtered_EC[as.Date(long_filtered_EC$Timestamp) != as.Date("2023-09-06"), ]


pdf("../Plots/Sensors/EC_all_growthhabit.pdf", height=5, width = 12)
p_EC <- ggplot(long_filtered_EC, aes(x=Timestamp, y=value, group = Sensor, colour = Sensor)) +
  geom_line() +
  geom_point(size = 1) + ylab("Bulk EC (mS/cm)") + xlab("Date") +
  theme_bw() + scale_color_brewer(palette = "Dark2", labels = c("Control determinate bush", "Drought determinate bush", "Drought indeterminate bush", "Drought indeterminate climbing")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_vline(xintercept = as.POSIXct("2023-07-27 00:00:00"), linetype="dotted", linewidth = 1) +
  geom_text(aes(x=as.POSIXct("2023-07-27 00:00:00"), label="\nIrrigation stopped", y=0.28), colour="black", angle=0, hjust = 0) +
  geom_vline(xintercept = as.POSIXct("2023-08-25 00:00:00"), linetype="dotted", linewidth = 1) +
  geom_text(aes(x= as.POSIXct("2023-08-25 00:00:00"), label="\nRecovery", y=0.28), colour="black", angle=0, hjust = 0) +
  scale_x_datetime(breaks = "3 days")

p_EC
dev.off()

#Violin plot EC
long_filtered_EC_more <- long_filtered_EC %>% filter(Timestamp >= "2023-07-27")

pdf("../Plots/Sensors/EC_all_violin_growthhabits_stats_irrigstopped.pdf", height=6, width = 8)
ps_EC <- ggstatsplot::ggbetweenstats(data = long_filtered_EC_more,
                                 x = Sensor,
                                 y = value,
                                 type = "parametric",
                                 pairwise.display = "significant",
                                 pairwise.annotation = "asterisk",
                                 p.adjust.method = "holm",
                                 results.subtitle = FALSE,
                                 centrality.point.args = list(size = 4, color = "darkred"),
                                 violin.args = list(width = 0.5, alpha = 0.2, na.rm = TRUE),
                                 ggsignif.args = list(textsize = 3, tip_length = 0.01, na.rm = TRUE),
                                 # mean.ci = TRUE,
                                 pairwise.comparisons = TRUE,
                                 bf.message = FALSE,
                                 xlab = "Sensors", 
                                 ylab = "Bulk EC (mS/cm)",
                                 digits = "signif",
                                 ggtheme = ggplot2::theme_bw(), 
                                 boxplot.args = list(width = 0) 
) +
  scale_x_discrete(labels = c(
    "Control_determinatebush" = "Control \nDeterminate Bush", 
    "Drought_determinatebush" = "Drought \nDeterminate Bush", 
    "Drought_indeterminatebush" = "Drought \nIndeterminate Bush", 
    "Drought_indeterminateclimbing" = "Drought \nIndeterminate Climbing"))

ps_EC
dev.off()

###Plotting Water content ----------------------------
#long_filtered <- long[c(-2988:-3030), ]
long_filtered <- long_WC[as.Date(long_WC$Timestamp) != as.Date("2023-09-08"), ]
long_filtered <- long_filtered[as.Date(long_filtered$Timestamp) != as.Date("2023-09-07"), ]
long_filtered <- long_filtered[as.Date(long_filtered$Timestamp) != as.Date("2023-09-06"), ]

pdf("../Plots/Sensors/water_content_all_growthhabit.pdf", height=5, width = 12)
p_WC <- ggplot(long_filtered, aes(x=Timestamp, y=value, group = Sensor, colour = Sensor)) +
  geom_line() +
  geom_point(size = 1) + ylab("Water content (m³/m³)") + xlab("Date") +
  theme_bw() + scale_color_brewer(palette = "Dark2", labels = c("Control determinate bush", "Drought determinate bush", "Drought indeterminate bush", "Drought indeterminate climbing")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_vline(xintercept = as.POSIXct("2023-07-27 00:00:00"), linetype="dotted", linewidth = 1) +
  geom_text(aes(x=as.POSIXct("2023-07-27 00:00:00"), label="\nIrrigation stopped", y=0.045), colour="black", angle=0, hjust = 0) +
  geom_vline(xintercept = as.POSIXct("2023-08-25 00:00:00"), linetype="dotted", linewidth = 1) +
  geom_text(aes(x= as.POSIXct("2023-08-25 00:00:00"), label="\nRecovery", y=0.045), colour="black", angle=0, hjust = 0) +
  scale_x_datetime(breaks = "3 days")
p_WC
dev.off()

#Violin plot water content
#Start from when irrigation was stopped for comparisons
long_filtered_more <- long_filtered %>% filter(Timestamp >= "2023-07-27")

pdf("../Plots/Sensors/water_content_all_violin_growthhabit_stats_irrigstopped.pdf", height=6, width = 8)
ps_WC <- ggstatsplot::ggbetweenstats(data = long_filtered_more,
                                 x = Sensor,
                                 y = value,
                                 type = "parametric",
                                 pairwise.display = "significant",
                                 pairwise.annotation = "asterisk",
                                 p.adjust.method = "holm",
                                 results.subtitle = FALSE,
                                 centrality.point.args = list(size = 4, color = "darkred"),
                                 violin.args = list(width = 0.5, alpha = 0.2, na.rm = TRUE),
                                 ggsignif.args = list(textsize = 3, tip_length = 0.01, na.rm = TRUE),
                                 # mean.ci = TRUE,
                                 pairwise.comparisons = TRUE,
                                 bf.message = FALSE,
                                 xlab = "Sensors", 
                                 ylab = "Water content (m³/m³)",
                                 digits = "signif",
                                 ggtheme = ggplot2::theme_bw(), 
                                 boxplot.args = list(width = 0) 
) +
  scale_x_discrete(labels = c(
    "Control_determinatebush" = "Control \nDeterminate Bush", 
    "Drought_determinatebush" = "Drought \nDeterminate Bush", 
    "Drought_indeterminatebush" = "Drought \nIndeterminate Bush", 
    "Drought_indeterminateclimbing" = "Drought \nIndeterminate Climbing"))


ps_WC
dev.off()

##Add vertical lines for removing water / recovery rain? 
#27/07/2023 water removed
# 25/08/2023 recovery

##### Put all of them in one plot  ----------
library(ggpubr)

pdf("../Plots/Sensors/allgrowthhabit_line.pdf", height=10, width = 12)
p1 <- ggarrange(p_MP, p_ST, p_EC, p_WC, 
                ncol = 2, 
                nrow = 2, 
                legend = "bottom",
                common.legend = T, 
                labels = "AUTO"
                )
p1
dev.off()


pdf("../Plots/Sensors/allgrowthhabit_violin_irrigstopped.pdf", height=11, width = 12)
p2 <- ggarrange(ps_MP, ps_ST, ps_EC, ps_WC, 
                ncol = 2, 
                nrow = 2, 
                legend = "none",
                labels = "AUTO"
)
p2
dev.off()
