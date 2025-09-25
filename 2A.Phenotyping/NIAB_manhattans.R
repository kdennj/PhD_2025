library(readxl) 
library(qqman)
library(ggplot2)
library(dplyr)
library(grDevices)
library(knitr)
library(ggpubr)
library(RColorBrewer)

setwd("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Bean project/23.03-Field-trial-NIAB/GAPIT/")

rm(list = ls())
gc()

set.seed(123) 

###GSW --------------------------------
#response_gsw_second (MLMM, BLINK, FarmCPU), response_gsw_diff_first_second (FarmCPU, BLINK),
#response_gsw_first (BLINK, FarmCPU)
#recovery_gsw_fifth (MLMM, BLINK, FarmCPU), recovery_gsw_fourth (FarmCPU, MLMM, BLINK),
#recovery_gsw_diff_fifth_fourth (MLMM, BLINK, FarmCPU)
#gsw_diff_first_fourth (BLINK, FarmCPU, MLMM), gsw_diff_first_third (BLINK, FarmCPU, MLMM)
#gsw_diff_first_fifth (MLMM, MLM, BLINK, FarmCPU)
#gsw_average (MLMM, FarmCPU, MLM, BLINK), 
#gsw_scaled_sum (MLMM, FarmCPU, MLM, BLINK), gsw_sum (MLMM, FarmCPU, MLM, BLINK), 

gsw_average_Res_MLMM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLMM.gsw_average.csv")
gsw_average_Res_FarmCPU <- read.csv("gsw/GAPIT.Association.GWAS_Results.FarmCPU.gsw_average.csv")
gsw_average_Res_MLM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLM.gsw_average.csv")
gsw_average_Res_BLINK <- read.csv("gsw/GAPIT.Association.GWAS_Results.BLINK.gsw_average.csv")
gsw_scaled_sum_Res_MLMM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLMM.gsw_scaled_sum.csv")
gsw_scaled_sum_Res_FarmCPU <- read.csv("gsw/GAPIT.Association.GWAS_Results.FarmCPU.gsw_scaled_sum.csv")
gsw_scaled_sum_Res_MLM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLM.gsw_scaled_sum.csv")
gsw_scaled_sum_Res_BLINK <- read.csv("gsw/GAPIT.Association.GWAS_Results.BLINK.gsw_scaled_sum.csv")
gsw_sum_Res_MLMM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLMM.gsw_sum.csv")
gsw_sum_Res_BLINK <- read.csv("gsw/GAPIT.Association.GWAS_Results.BLINK.gsw_sum.csv")
gsw_sum_Res_MLM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLM.gsw_sum.csv")
gsw_sum_Res_FarmCPU <- read.csv("gsw/GAPIT.Association.GWAS_Results.FarmCPU.gsw_sum.csv")

rm(list = ls(pattern = "response"))

gsw_diff_first_fourth_Res_BLINK <- read.csv("gsw/GAPIT.Association.GWAS_Results.BLINK.gsw_diff_first_fourth.csv")
gsw_diff_first_fourth_Res_FarmCPU <- read.csv("gsw/GAPIT.Association.GWAS_Results.FarmCPU.gsw_diff_first_fourth.csv")
gsw_diff_first_fourth_Res_MLMM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLMM.gsw_diff_first_fourth.csv")
gsw_diff_first_third_Res_BLINK <- read.csv("gsw/GAPIT.Association.GWAS_Results.BLINK.gsw_diff_first_third.csv")
gsw_diff_first_third_Res_FarmCPU <- read.csv("gsw/GAPIT.Association.GWAS_Results.FarmCPU.gsw_diff_first_third.csv")
gsw_diff_first_third_Res_MLMM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLMM.gsw_diff_first_third.csv")
gsw_diff_first_fifth_Res_BLINK <- read.csv("gsw/GAPIT.Association.GWAS_Results.BLINK.gsw_diff_first_fifth.csv")
gsw_diff_first_fifth_Res_FarmCPU <- read.csv("gsw/GAPIT.Association.GWAS_Results.FarmCPU.gsw_diff_first_fifth.csv")
gsw_diff_first_fifth_Res_MLMM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLMM.gsw_diff_first_fifth.csv")
#Bad QQ plots
gsw_diff_first_fourth_Res_MLM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLM.gsw_diff_first_fourth.csv")
gsw_diff_first_third_Res_MLM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLM.gsw_diff_first_third.csv")
gsw_diff_first_fifth_Res_MLM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLM.gsw_diff_first_fifth.csv")


rm(list = ls(pattern = "leaf"))


### Transp ------------------------------
#response_transp_first (BLINK, FarmCPU)
#response_transp_second (BLINK, FarmCPU), recovery_transp_fourth (BLINK, MLMM, FarmCPU)
#Transp_average (BLINK, MLMM, FarmCPU), Transp_diff_first_fifth (FarmCPU)
#Transp_diff_first_fourth (BLINK, MLMM, FarmCPU), Transp_diff_first_third (BLINK, MLMM, FarmCPU)
#Transp_scaled_sum (BLINK, MLMM, FarmCPU), Transp_sum (BLINK, MLMM, FarmCPU)

rm(list = ls(pattern = "Transp"))

Transp_diff_first_fifth_Res_FarmCPU <- read.csv("Transp/GAPIT.Association.GWAS_Results.FarmCPU.Transp_diff_first_fifth.csv")
Transp_diff_first_fourth_Res_FarmCPU <- read.csv("Transp/GAPIT.Association.GWAS_Results.FarmCPU.Transp_diff_first_fourth.csv")
Transp_diff_first_fourth_Res_BLINK <- read.csv("Transp/GAPIT.Association.GWAS_Results.BLINK.Transp_diff_first_fourth.csv")
Transp_diff_first_fourth_Res_MLMM <- read.csv("Transp/GAPIT.Association.GWAS_Results.MLMM.Transp_diff_first_fourth.csv")
Transp_diff_first_third_Res_MLMM <- read.csv("Transp/GAPIT.Association.GWAS_Results.MLMM.Transp_diff_first_third.csv")
Transp_diff_first_third_Res_BLINK <- read.csv("Transp/GAPIT.Association.GWAS_Results.BLINK.Transp_diff_first_third.csv")
Transp_diff_first_third_Res_FarmCPU <- read.csv("Transp/GAPIT.Association.GWAS_Results.FarmCPU.Transp_diff_first_third.csv")
#Bad QQ plots
Transp_diff_first_fifth_Res_MLM <- read.csv("Transp/GAPIT.Association.GWAS_Results.MLM.Transp_diff_first_fifth.csv")
Transp_diff_first_fifth_Res_MLMM <- read.csv("Transp/GAPIT.Association.GWAS_Results.MLMM.Transp_diff_first_fifth.csv")
Transp_diff_first_fifth_Res_BLINK <- read.csv("Transp/GAPIT.Association.GWAS_Results.BLINK.Transp_diff_first_fifth.csv")
Transp_diff_first_fourth_Res_MLM <- read.csv("Transp/GAPIT.Association.GWAS_Results.MLM.Transp_diff_first_fourth.csv")
Transp_diff_first_third_Res_MLM <- read.csv("Transp/GAPIT.Association.GWAS_Results.MLM.Transp_diff_first_third.csv")


Transp_average_Res_BLINK <- read.csv("Transp/GAPIT.Association.GWAS_Results.BLINK.Transp_average.csv")
Transp_average_Res_MLMM <- read.csv("Transp/GAPIT.Association.GWAS_Results.MLMM.Transp_average.csv")
Transp_average_Res_FarmCPU <- read.csv("Transp/GAPIT.Association.GWAS_Results.FarmCPU.Transp_average.csv")
Transp_scaled_sum_Res_BLINK <- read.csv("Transp/GAPIT.Association.GWAS_Results.BLINK.Transp_scaled_sum.csv")
Transp_scaled_sum_Res_FarmCPU <- read.csv("Transp/GAPIT.Association.GWAS_Results.FarmCPU.Transp_scaled_sum.csv")
Transp_scaled_sum_Res_MLMM <- read.csv("Transp/GAPIT.Association.GWAS_Results.MLMM.Transp_scaled_sum.csv")
Transp_sum_Res_MLMM <- read.csv("Transp/GAPIT.Association.GWAS_Results.MLMM.Transp_sum.csv")
Transp_sum_Res_FarmCPU <- read.csv("Transp/GAPIT.Association.GWAS_Results.FarmCPU.Transp_sum.csv")
Transp_sum_Res_BLINK <- read.csv("Transp/GAPIT.Association.GWAS_Results.BLINK.Transp_sum.csv")
#Bad QQ plots
Transp_average_Res_MLM <- read.csv("Transp/GAPIT.Association.GWAS_Results.MLM.Transp_average.csv")
Transp_scaled_sum_Res_MLM <- read.csv("Transp/GAPIT.Association.GWAS_Results.MLM.Transp_scaled_sum.csv")
Transp_sum_Res_MLM <- read.csv("Transp/GAPIT.Association.GWAS_Results.MLM.Transp_sum.csv")


###BBCH ------------------------------
#BBCH_average (BLINK, FarmCPU), BBCH_diff_fifth_first (BLINK, FarmCPU),
#BBCH_diff_fourth_first (BLINK, FarmCPU), BBCH_diff_third_first (BLINK),
#BBCH_sum (BLINK, FarmCPU)

rm(list = ls(pattern = "ETR"))

BBCH_average_Res_BLINK <- read.csv("BBCH/GAPIT.Association.GWAS_Results.BLINK.BBCH_average.csv")
BBCH_average_Res_FarmCPU <- read.csv("BBCH/GAPIT.Association.GWAS_Results.FarmCPU.BBCH_average.csv")
BBCH_diff_fifth_first_Res_BLINK <- read.csv("BBCH/GAPIT.Association.GWAS_Results.BLINK.BBCH_diff_fifth_first.csv")
BBCH_diff_fifth_first_Res_FarmCPU <- read.csv("BBCH/GAPIT.Association.GWAS_Results.FarmCPU.BBCH_diff_fifth_first.csv")
BBCH_diff_fourth_first_Res_BLINK <- read.csv("BBCH/GAPIT.Association.GWAS_Results.BLINK.BBCH_diff_fourth_first.csv")
BBCH_diff_fourth_first_Res_FarmCPU <- read.csv("BBCH/GAPIT.Association.GWAS_Results.FarmCPU.BBCH_diff_fourth_first.csv")
BBCH_diff_third_first_Res_BLINK <- read.csv("BBCH/GAPIT.Association.GWAS_Results.BLINK.BBCH_diff_third_first.csv")
BBCH_sum_Res_BLINK <- read.csv("BBCH/GAPIT.Association.GWAS_Results.BLINK.BBCH_sum.csv")
BBCH_sum_Res_FarmCPU <- read.csv("BBCH/GAPIT.Association.GWAS_Results.FarmCPU.BBCH_sum.csv")
#Bad QQ plots
BBCH_average_Res_MLM <- read.csv("BBCH/GAPIT.Association.GWAS_Results.MLM.BBCH_average.csv")
BBCH_average_Res_MLMM <- read.csv("BBCH/GAPIT.Association.GWAS_Results.MLMM.BBCH_average.csv")
BBCH_diff_fifth_first_Res_MLM <- read.csv("BBCH/GAPIT.Association.GWAS_Results.MLM.BBCH_diff_fifth_first.csv")
BBCH_diff_fifth_first_Res_MLMM <- read.csv("BBCH/GAPIT.Association.GWAS_Results.MLMM.BBCH_diff_fifth_first.csv")
BBCH_diff_fourth_first_Res_MLM <- read.csv("BBCH/GAPIT.Association.GWAS_Results.MLM.BBCH_diff_fourth_first.csv")
BBCH_diff_fourth_first_Res_MLMM <- read.csv("BBCH/GAPIT.Association.GWAS_Results.MLMM.BBCH_diff_fourth_first.csv")
BBCH_diff_third_first_Res_FarmCPU <- read.csv("BBCH/GAPIT.Association.GWAS_Results.FarmCPU.BBCH_diff_third_first.csv")
BBCH_diff_third_first_Res_MLM <- read.csv("BBCH/GAPIT.Association.GWAS_Results.MLM.BBCH_diff_third_first.csv")
BBCH_diff_third_first_Res_MLMM <- read.csv("BBCH/GAPIT.Association.GWAS_Results.MLMM.BBCH_diff_third_first.csv")
BBCH_sum_Res_MLM <- read.csv("BBCH/GAPIT.Association.GWAS_Results.MLM.BBCH_sum.csv")
BBCH_sum_Res_MLMM <- read.csv("BBCH/GAPIT.Association.GWAS_Results.MLMM.BBCH_sum.csv")



###ETR ------------------------------
#recovery_ETR_fifth (BLINK, MLMM, FarmCPU), 
#ETR_diff_first_fifth (MLMM, BLINK, FarmCPU, MLM)

rm(list = ls(pattern = "weight"))

ETR_diff_first_fifth_Res_MLMM <- read.csv("ETR/GAPIT.Association.GWAS_Results.MLMM.ETR_diff_first_fifth.csv")
ETR_diff_first_fifth_Res_BLINK <- read.csv("ETR/GAPIT.Association.GWAS_Results.BLINK.ETR_diff_first_fifth.csv")
ETR_diff_first_fifth_Res_FarmCPU <- read.csv("ETR/GAPIT.Association.GWAS_Results.FarmCPU.ETR_diff_first_fifth.csv")
ETR_diff_first_fifth_Res_MLM <- read.csv("ETR/GAPIT.Association.GWAS_Results.MLM.ETR_diff_first_fifth.csv")



###Harvesting ------------------------------
#No_pods (FarmCPU, BLINK), 
#Partitioning_weight (FarmCPU, BLINK), weight_per_pod (FarmCPU), 
#Weight_pods (FarmCPU)

rm(list = ls(pattern = "BBCH"))

No_pods_Res_BLINK <- read.csv("Harvesting/GAPIT.Association.GWAS_Results.BLINK.No_pods.csv")
No_pods_Res_FarmCPU <- read.csv("Harvesting/GAPIT.Association.GWAS_Results.FarmCPU.No_pods.csv")
Partitioning_weight_Res_BLINK <- read.csv("Harvesting/GAPIT.Association.GWAS_Results.BLINK.Partitioning_weight.csv")
Partitioning_weight_Res_FarmCPU <- read.csv("Harvesting/GAPIT.Association.GWAS_Results.FarmCPU.Partitioning_weight.csv")
Weight_per_pod_Res_FarmCPU <- read.csv("Harvesting/GAPIT.Association.GWAS_Results.FarmCPU.weight_per_pod.csv")
Weight_pods_Res_FarmCPU <- read.csv("Harvesting/GAPIT.Association.GWAS_Results.FarmCPU.weight_pods.csv")
#Bad QQ plots
No_pods_Res_MLM <- read.csv("Harvesting/GAPIT.Association.GWAS_Results.MLM.No_pods.csv")
No_pods_Res_MLMM <- read.csv("Harvesting/GAPIT.Association.GWAS_Results.MLMM.No_pods.csv")
Partitioning_weight_Res_MLM <- read.csv("Harvesting/GAPIT.Association.GWAS_Results.MLM.Partitioning_weight.csv")
#Partitioning_weight_Res_MLMM <- read.csv("Harvesting/GAPIT.Association.GWAS_Results.MLMM.Partitioning_weight.csv")
Weight_per_pod_Res_BLINK <- read.csv("Harvesting/GAPIT.Association.GWAS_Results.BLINK.weight_per_pod.csv")
Weight_per_pod_Res_MLM <- read.csv("Harvesting/GAPIT.Association.GWAS_Results.MLM.weight_per_pod.csv")
#Weight_per_pod_Res_MLMM <- read.csv("Harvesting/GAPIT.Association.GWAS_Results.MLMM.weight_per_pod.csv")
Weight_pods_Res_BLINK <- read.csv("Harvesting/GAPIT.Association.GWAS_Results.BLINK.weight_pods.csv")
Weight_pods_Res_MLM <- read.csv("Harvesting/GAPIT.Association.GWAS_Results.MLM.weight_pods.csv")
#Weight_pods_Res_MLMM <- read.csv("Harvesting/GAPIT.Association.GWAS_Results.MLMM.weight_pods.csv")



###VPDleaf ------------------------------
#response_VPDleaf_first (BLINK), VPDleaf_diff_first_fourth (FarmCPU, BLINK, MLMM), 
#VPDleaf_diff_first_third (MLMM)

rm(list = ls(pattern = "ETR"))

VPDleaf_diff_first_fourth_Res_FarmCPU <- read.csv("VPDleaf/GAPIT.Association.GWAS_Results.FarmCPU.VPDleaf_diff_first_fourth.csv")
VPDleaf_diff_first_fourth_Res_BLINK <- read.csv("VPDleaf/GAPIT.Association.GWAS_Results.BLINK.VPDleaf_diff_first_fourth.csv")
VPDleaf_diff_first_fourth_Res_MLMM <- read.csv("VPDleaf/GAPIT.Association.GWAS_Results.MLMM.VPDleaf_diff_first_fourth.csv")
VPDleaf_diff_first_third_Res_MLMM <- read.csv("VPDleaf/GAPIT.Association.GWAS_Results.MLMM.VPDleaf_diff_first_third.csv")
#Bad QQ plots
VPDleaf_diff_first_fourth_Res_MLM <- read.csv("VPDleaf/GAPIT.Association.GWAS_Results.MLM.VPDleaf_diff_first_fourth.csv")
VPDleaf_diff_first_third_Res_MLM <- read.csv("VPDleaf/GAPIT.Association.GWAS_Results.MLM.VPDleaf_diff_first_third.csv")
VPDleaf_diff_first_third_Res_FarmCPU <- read.csv("VPDleaf/GAPIT.Association.GWAS_Results.FarmCPU.VPDleaf_diff_first_third.csv")
VPDleaf_diff_first_third_Res_BLINK <- read.csv("VPDleaf/GAPIT.Association.GWAS_Results.BLINK.VPDleaf_diff_first_third.csv")


###Tleaf ------------------------------
#Tleaf_diff_first_third (MLMM), Tleaf_diff_first_fourth (FarmCPU, MLMM)

Tleaf_first_third_Res_MLMM <- read.csv("Tleaf/GAPIT.Association.GWAS_Results.MLMM.Tleaf_diff_first_third.csv")
Tleaf_first_fourth_Res_MLMM <- read.csv("Tleaf/GAPIT.Association.GWAS_Results.MLMM.Tleaf_diff_first_fourth.csv")
Tleaf_first_fourth_Res_FarmCPU <- read.csv("Tleaf/GAPIT.Association.GWAS_Results.FarmCPU.Tleaf_diff_first_fourth.csv")
#Bad QQ plots
Tleaf_first_third_Res_MLM <- read.csv("Tleaf/GAPIT.Association.GWAS_Results.MLM.Tleaf_diff_first_third.csv")
Tleaf_first_third_Res_BLINK <- read.csv("Tleaf/GAPIT.Association.GWAS_Results.BLINK.Tleaf_diff_first_third.csv")
Tleaf_first_third_Res_FarmCPU <- read.csv("Tleaf/GAPIT.Association.GWAS_Results.FarmCPU.Tleaf_diff_first_third.csv")
Tleaf_first_fourth_Res_MLM <- read.csv("Tleaf/GAPIT.Association.GWAS_Results.MLM.Tleaf_diff_first_fourth.csv")
Tleaf_first_fourth_Res_BLINK <- read.csv("Tleaf/GAPIT.Association.GWAS_Results.BLINK.Tleaf_diff_first_fourth.csv")


###Response ------------------------------

VPDleaf_response_first_Res_BLINK <- read.csv("VPDleaf/GAPIT.Association.GWAS_Results.BLINK.response_VPDleaf_first.csv")
response_transp_first_Res_BLINK <- read.csv("Transp/GAPIT.Association.GWAS_Results.BLINK.response_transp_first.csv")
response_transp_first_Res_FarmCPU <- read.csv("Transp/GAPIT.Association.GWAS_Results.FarmCPU.response_transp_first.csv")
response_transp_second_Res_BLINK <- read.csv("Transp/GAPIT.Association.GWAS_Results.BLINK.response_transp_second.csv")
response_transp_second_Res_FarmCPU <- read.csv("Transp/GAPIT.Association.GWAS_Results.FarmCPU.response_transp_second.csv")
response_gsw_second_Res_MLMM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLMM.response_gsw_second.csv")
response_gsw_second_Res_BLINK <- read.csv("gsw/GAPIT.Association.GWAS_Results.BLINK.response_gsw_second.csv")
response_gsw_second_Res_FarmCPU <- read.csv("gsw/GAPIT.Association.GWAS_Results.FarmCPU.response_gsw_second.csv")
response_gsw_diff_first_second_Res_BLINK <- read.csv("gsw/GAPIT.Association.GWAS_Results.BLINK.response_gsw_diff_first_second.csv")
response_gsw_diff_first_second_Res_FarmCPU <- read.csv("gsw/GAPIT.Association.GWAS_Results.FarmCPU.response_gsw_diff_first_second.csv")
response_gsw_first_Res_BLINK <- read.csv("gsw/GAPIT.Association.GWAS_Results.BLINK.response_gsw_first.csv")
response_gsw_first_Res_FarmCPU <- read.csv("gsw/GAPIT.Association.GWAS_Results.FarmCPU.response_gsw_first.csv")
#Bad QQ plots
VPDleaf_response_first_Res_FarmCPU <- read.csv("VPDleaf/GAPIT.Association.GWAS_Results.FarmCPU.response_VPDleaf_first.csv")
VPDleaf_response_first_Res_MLM <- read.csv("VPDleaf/GAPIT.Association.GWAS_Results.MLM.response_VPDleaf_first.csv")
VPDleaf_response_first_Res_MLMM <- read.csv("VPDleaf/GAPIT.Association.GWAS_Results.MLMM.response_VPDleaf_first.csv")
response_transp_first_Res_MLM <- read.csv("Transp/GAPIT.Association.GWAS_Results.MLM.response_transp_first.csv")
response_transp_first_Res_MLMM <- read.csv("Transp/GAPIT.Association.GWAS_Results.MLMM.response_transp_first.csv")
response_transp_second_Res_MLM <- read.csv("Transp/GAPIT.Association.GWAS_Results.MLM.response_transp_second.csv")
response_transp_second_Res_MLMM <- read.csv("Transp/GAPIT.Association.GWAS_Results.MLMM.response_transp_second.csv")
response_gsw_second_Res_MLM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLM.response_gsw_second.csv")
response_gsw_diff_first_second_Res_MLM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLM.response_gsw_diff_first_second.csv")
response_gsw_diff_first_second_Res_MLMM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLMM.response_gsw_diff_first_second.csv")
response_gsw_first_Res_MLM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLM.response_gsw_first.csv")
response_gsw_first_Res_MLMM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLMM.response_gsw_first.csv")


###Recovery ------------------------------

recovery_ETR_fifth_Res_BLINK <- read.csv("ETR/GAPIT.Association.GWAS_Results.BLINK.recovery_ETR_fifth.csv")
recovery_ETR_fifth_Res_MLMM <- read.csv("ETR/GAPIT.Association.GWAS_Results.MLMM.recovery_ETR_fifth.csv")
recovery_ETR_fifth_Res_FarmCPU <- read.csv("ETR/GAPIT.Association.GWAS_Results.FarmCPU.recovery_ETR_fifth.csv")
recovery_transp_fourth_Res_BLINK <- read.csv("Transp/GAPIT.Association.GWAS_Results.BLINK.recovery_transp_fourth.csv")
recovery_transp_fourth_Res_MLMM <- read.csv("Transp/GAPIT.Association.GWAS_Results.MLMM.recovery_transp_fourth.csv")
recovery_transp_fourth_Res_FarmCPU <- read.csv("Transp/GAPIT.Association.GWAS_Results.FarmCPU.recovery_transp_fourth.csv")
recovery_gsw_fifth_Res_MLMM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLMM.recovery_gsw_fifth.csv")
recovery_gsw_fifth_Res_BLINK <- read.csv("gsw/GAPIT.Association.GWAS_Results.BLINK.recovery_gsw_fifth.csv")
recovery_gsw_fifth_Res_FarmCPU <- read.csv("gsw/GAPIT.Association.GWAS_Results.FarmCPU.recovery_gsw_fifth.csv")
recovery_gsw_fourth_Res_MLMM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLMM.recovery_gsw_fourth.csv")
recovery_gsw_fourth_Res_BLINK <- read.csv("gsw/GAPIT.Association.GWAS_Results.BLINK.recovery_gsw_fourth.csv")
recovery_gsw_fourth_Res_FarmCPU <- read.csv("gsw/GAPIT.Association.GWAS_Results.FarmCPU.recovery_gsw_fourth.csv")
recovery_gsw_diff_fifth_fourth_Res_MLMM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLMM.recovery_gsw_diff_fifth_fourth.csv")
recovery_gsw_diff_fifth_fourth_Res_BLINK <- read.csv("gsw/GAPIT.Association.GWAS_Results.BLINK.recovery_gsw_diff_fifth_fourth.csv")
recovery_gsw_diff_fifth_fourth_Res_FarmCPU <- read.csv("gsw/GAPIT.Association.GWAS_Results.FarmCPU.recovery_gsw_diff_fifth_fourth.csv")
#Bad QQ plots
recovery_ETR_fifth_Res_MLM <- read.csv("ETR/GAPIT.Association.GWAS_Results.MLM.recovery_ETR_fifth.csv")
recovery_transp_fourth_Res_MLM <- read.csv("Transp/GAPIT.Association.GWAS_Results.MLM.recovery_transp_fourth.csv")
recovery_gsw_fifth_Res_MLM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLM.recovery_gsw_fifth.csv")
recovery_gsw_fourth_Res_MLM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLM.recovery_gsw_fourth.csv")
recovery_gsw_diff_fifth_fourth_Res_MLM <- read.csv("gsw/GAPIT.Association.GWAS_Results.MLM.recovery_gsw_diff_fifth_fourth.csv")


##### Add Model column ------------------------------
# df_list <- list(
#   Tleaf_first_fourth_Res_FarmCPU = Tleaf_first_fourth_Res_FarmCPU, 
#   Tleaf_first_fourth_Res_MLMM = Tleaf_first_fourth_Res_MLMM, 
#   Tleaf_first_third_Res_MLMM = Tleaf_first_third_Res_MLMM
# )
# 
# for (df_name in names(df_list)) {
#   # Extract last part of the name using regex
#   model_name <- sub(".*_(\\w+)$", "\\1", df_name)
#   
#   # Add model column
#   df_list[[df_name]]$model <- model_name
# }

# Combine the dataframes
# gwasResults <- do.call(rbind, df_list)

####### Basic Manhattan with ggplot ----------------------------------------------------------------------------
results <- df_list$Tleaf_first_fourth_Res_FarmCPU
results <- Tleaf_first_third_Res_MLMM

don <- results %>% 
  
  # Compute chromosome size
  group_by(Chr) %>% 
  summarise(chr_len=max(Pos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(results, ., by=c("Chr"="Chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chr, Pos) %>%
  mutate( BPcum=Pos+tot)

axisdf = don %>%
  group_by(Chr) %>%
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

sig <- -log10(0.05 / nrow(results))

#Plot with ggplot2 
par(xaxs = "i", yaxs = "i")

p <- ggplot(don, aes(x=BPcum, y=-log10(P.value))) +
  
  # Show all points
  geom_point( aes(color=as.factor(Chr)), alpha=0.8, size=2) +
  #scale_color_manual(values = rep(c("#ff6666", "#ffbd55", "#9de24f", "#87cefa", "#cf9dff"), 11 )) +
  scale_color_manual(values = rep(c("black", "darkgrey"), 11 )) +
  
  #Bonferonni cutoff
  geom_hline(
    yintercept = sig, color = "forestgreen") +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$Chr, breaks= axisdf$center ) +
  #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_classic() +
  theme( 
    axis.title.x=element_blank(),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank() 
  )
#p

#pdf(file ="det_A_gwas_F_manhattan.pdf", width = 10, height = 4)
#plot(p)
#dev.off()

jpeg("Tleaf/Tleaf_first_fourth_Res_FarmCPU.jpg", width = 1500, height = 600)
plot(p)
dev.off()

##### Combine manhattan plots --------------------------------------------------

# VPDleaf, Tleaf and ETR
df_list <- list(
  VPDleaf_diff_first_fourth_Res_BLINK = VPDleaf_diff_first_fourth_Res_BLINK,
  VPDleaf_diff_first_fourth_Res_FarmCPU = VPDleaf_diff_first_fourth_Res_FarmCPU,
  VPDleaf_diff_first_fourth_Res_MLMM = VPDleaf_diff_first_fourth_Res_MLMM,
  VPDleaf_diff_first_fourth_Res_MLM = VPDleaf_diff_first_fourth_Res_MLM,
  VPDleaf_diff_first_third_Res_MLMM = VPDleaf_diff_first_third_Res_MLMM,
  VPDleaf_diff_first_third_Res_MLM = VPDleaf_diff_first_third_Res_MLM,
  VPDleaf_diff_first_third_Res_FarmCPU = VPDleaf_diff_first_third_Res_FarmCPU,
  VPDleaf_diff_first_third_Res_BLINK = VPDleaf_diff_first_third_Res_BLINK,
  Tleaf_diff_first_fourth_Res_FarmCPU = Tleaf_first_fourth_Res_FarmCPU,
  Tleaf_diff_first_fourth_Res_MLMM = Tleaf_first_fourth_Res_MLMM,
  Tleaf_diff_first_fourth_Res_MLM = Tleaf_first_fourth_Res_MLM,
  Tleaf_diff_first_fourth_Res_BLINK = Tleaf_first_fourth_Res_BLINK,
  Tleaf_diff_first_third_Res_MLMM = Tleaf_first_third_Res_MLMM,
  Tleaf_diff_first_third_Res_MLM = Tleaf_first_third_Res_MLM,
  Tleaf_diff_first_third_Res_BLINK = Tleaf_first_third_Res_BLINK,
  Tleaf_diff_first_third_Res_FarmCPU = Tleaf_first_third_Res_FarmCPU,
  ETR_diff_first_fifth_Res_BLINK = ETR_diff_first_fifth_Res_BLINK,
  ETR_diff_first_fifth_Res_FarmCPU = ETR_diff_first_fifth_Res_FarmCPU,
  ETR_diff_first_fifth_Res_MLM = ETR_diff_first_fifth_Res_MLM,
  ETR_diff_first_fifth_Res_MLMM = ETR_diff_first_fifth_Res_MLMM
)

#Harvesting
df_list <- list(
  No_pods_Res_BLINK = No_pods_Res_BLINK,
  No_pods_Res_FarmCPU = No_pods_Res_FarmCPU,
  Partitioning_weight_Res_BLINK = Partitioning_weight_Res_BLINK,
  Partitioning_weight_Res_FarmCPU = Partitioning_weight_Res_FarmCPU,
  Weight_per_pod_Res_FarmCPU = Weight_per_pod_Res_FarmCPU,
  Weight_pods_Res_FarmCPU = Weight_pods_Res_FarmCPU,
  
  No_pods_Res_MLM = No_pods_Res_MLM,
  No_pods_Res_MLMM = No_pods_Res_MLMM,
  Partitioning_weight_Res_MLM = Partitioning_weight_Res_MLM,
  Weight_per_pod_Res_BLINK = Weight_per_pod_Res_BLINK,
  Weight_per_pod_Res_MLM = Weight_per_pod_Res_MLM,
  Weight_pods_Res_BLINK = Weight_pods_Res_BLINK,
  Weight_pods_Res_MLM = Weight_pods_Res_MLM
)

#BBCH
df_list <- list(
  BBCH_average_Res_BLINK = BBCH_average_Res_BLINK,
  BBCH_average_Res_FarmCPU = BBCH_average_Res_FarmCPU,
  BBCH_diff_fifth_first_Res_BLINK = BBCH_diff_fifth_first_Res_BLINK,
  BBCH_diff_fifth_first_Res_FarmCPU = BBCH_diff_fifth_first_Res_FarmCPU,
  BBCH_diff_fourth_first_Res_BLINK = BBCH_diff_fourth_first_Res_BLINK,
  BBCH_diff_fourth_first_Res_FarmCPU = BBCH_diff_fourth_first_Res_FarmCPU,
  BBCH_diff_third_first_Res_BLINK = BBCH_diff_third_first_Res_BLINK,
  BBCH_sum_Res_BLINK = BBCH_sum_Res_BLINK,
  BBCH_sum_Res_FarmCPU = BBCH_sum_Res_FarmCPU,
  
  BBCH_average_Res_MLM = BBCH_average_Res_MLM,
  BBCH_average_Res_MLMM = BBCH_average_Res_MLMM,
  BBCH_diff_fifth_first_Res_MLMM = BBCH_diff_fifth_first_Res_MLMM,
  BBCH_diff_fifth_first_Res_MLM = BBCH_diff_fifth_first_Res_MLM,
  BBCH_diff_fourth_first_Res_MLMM = BBCH_diff_fourth_first_Res_MLMM ,
  BBCH_diff_fourth_first_Res_MLM = BBCH_diff_fourth_first_Res_MLM,
  BBCH_diff_third_first_Res_FarmCPU = BBCH_diff_third_first_Res_FarmCPU,
  BBCH_diff_third_first_Res_MLM = BBCH_diff_third_first_Res_MLM,
  BBCH_diff_third_first_Res_MLMM = BBCH_diff_third_first_Res_MLMM,
  BBCH_sum_Res_MLM = BBCH_sum_Res_MLM,
  BBCH_sum_Res_MLMM = BBCH_sum_Res_MLMM
)

#Transp_response_recovery
# df_list <- list(
#   recovery_transp_fourth_Res_BLINK = recovery_transp_fourth_Res_BLINK,
#   recovery_transp_fourth_Res_FarmCPU = recovery_transp_fourth_Res_FarmCPU,
#   recovery_transp_fourth_Res_MLMM = recovery_transp_fourth_Res_MLMM,
#   response_transp_first_Res_BLINK = response_transp_first_Res_BLINK,
#   response_transp_first_Res_FarmCPU = response_transp_first_Res_FarmCPU,
#   response_transp_second_Res_BLINK = response_transp_second_Res_BLINK,
#   response_transp_second_Res_FarmCPU = response_transp_second_Res_FarmCPU
# )

#Transp diff
df_list <- list(
  Transp_diff_first_fifth_Res_FarmCPU = Transp_diff_first_fifth_Res_FarmCPU,
  Transp_diff_first_fifth_Res_MLMM = Transp_diff_first_fifth_Res_MLMM,
  Transp_diff_first_fifth_Res_MLM = Transp_diff_first_fifth_Res_MLM,
  Transp_diff_first_fifth_Res_BLINK = Transp_diff_first_fifth_Res_BLINK,
  Transp_diff_first_fourth_Res_BLINK = Transp_diff_first_fourth_Res_BLINK,
  Transp_diff_first_fourth_Res_FarmCPU = Transp_diff_first_fourth_Res_FarmCPU,
  Transp_diff_first_fourth_Res_MLMM = Transp_diff_first_fourth_Res_MLMM,
  Transp_diff_first_fourth_Res_MLM = Transp_diff_first_fourth_Res_MLM,
  Transp_diff_first_third_Res_BLINK = Transp_diff_first_third_Res_BLINK,
  Transp_diff_first_third_Res_FarmCPU = Transp_diff_first_third_Res_FarmCPU,
  Transp_diff_first_third_Res_MLMM = Transp_diff_first_third_Res_MLMM,
  Transp_diff_first_third_Res_MLM = Transp_diff_first_third_Res_MLM
)

#Transp
df_list <- list(
  Transp_average_Res_BLINK = Transp_average_Res_BLINK,
  Transp_average_Res_FarmCPU = Transp_average_Res_FarmCPU,
  Transp_average_Res_MLMM = Transp_average_Res_MLMM,
  Transp_average_Res_MLM = Transp_average_Res_MLM,
  Transp_scaled_sum_Res_BLINK = Transp_scaled_sum_Res_BLINK,
  Transp_scaled_sum_Res_FarmCPU = Transp_scaled_sum_Res_FarmCPU,
  Transp_scaled_sum_Res_MLMM = Transp_scaled_sum_Res_MLMM,
  Transp_scaled_sum_Res_MLM = Transp_scaled_sum_Res_MLM,
  Transp_sum_Res_BLINK = Transp_sum_Res_BLINK,
  Transp_sum_Res_FarmCPU = Transp_sum_Res_FarmCPU,
  Transp_sum_Res_MLMM = Transp_scaled_sum_Res_MLMM,
  Transp_sum_Res_MLM = Transp_scaled_sum_Res_MLM
)

#gsw response
# df_list <- list(
#   response_gsw_diff_first_second_Res_BLINK = response_gsw_diff_first_second_Res_BLINK,
#   response_gsw_diff_first_second_Res_FarmCPU = response_gsw_diff_first_second_Res_FarmCPU,
#   response_gsw_first_Res_BLINK = response_gsw_first_Res_BLINK,
#   response_gsw_first_Res_FarmCPU = response_gsw_first_Res_FarmCPU,
#   response_gsw_second_Res_BLINK = response_gsw_second_Res_BLINK,
#   response_gsw_second_Res_FarmCPU = response_gsw_second_Res_FarmCPU,
#   response_gsw_second_Res_MLMM = response_gsw_second_Res_MLMM
# )

#gsw_recovery
# df_list <- list(
#   recovery_gsw_diff_fifth_fourth_Res_BLINK = recovery_gsw_diff_fifth_fourth_Res_BLINK,
#   recovery_gsw_diff_fifth_fourth_Res_FarmCPU = recovery_gsw_diff_fifth_fourth_Res_FarmCPU,
#   recovery_gsw_diff_fifth_fourth_Res_MLMM = recovery_gsw_diff_fifth_fourth_Res_MLMM,
#   recovery_gsw_fifth_Res_BLINK = recovery_gsw_fifth_Res_BLINK,
#   recovery_gsw_fifth_Res_FarmCPU = recovery_gsw_fifth_Res_FarmCPU,
#   recovery_gsw_fifth_Res_MLMM = recovery_gsw_fifth_Res_MLMM,
#   recovery_gsw_fourth_Res_BLINK = recovery_gsw_fourth_Res_BLINK,
#   recovery_gsw_fourth_Res_FarmCPU = recovery_gsw_fourth_Res_FarmCPU,
#   recovery_gsw_fourth_Res_MLMM = recovery_gsw_fourth_Res_MLMM
# )

#gsw diff - with QQ plot
df_list <- list(
  gsw_diff_first_fifth_Res_BLINK = gsw_diff_first_fifth_Res_BLINK,
  gsw_diff_first_fifth_Res_FarmCPU = gsw_diff_first_fifth_Res_FarmCPU,
  gsw_diff_first_fifth_Res_MLMM = gsw_diff_first_fifth_Res_MLMM,
  gsw_diff_first_fifth_Res_MLM = gsw_diff_first_fifth_Res_MLM,
  gsw_diff_first_fourth_Res_BLINK = gsw_diff_first_fourth_Res_BLINK,
  gsw_diff_first_fourth_Res_FarmCPU = gsw_diff_first_fourth_Res_FarmCPU,
  gsw_diff_first_fourth_Res_MLMM = gsw_diff_first_fourth_Res_MLMM,
  gsw_diff_first_fourth_Res_MLM = gsw_diff_first_fourth_Res_MLM,
  gsw_diff_first_third_Res_BLINK = gsw_diff_first_third_Res_BLINK,
  gsw_diff_first_third_Res_FarmCPU = gsw_diff_first_third_Res_FarmCPU,
  gsw_diff_first_third_Res_MLMM = gsw_diff_first_third_Res_MLMM,
  gsw_diff_first_third_Res_MLM = gsw_diff_first_third_Res_MLM
)

#gsw
df_list <- list(
  gsw_average_Res_BLINK = gsw_average_Res_BLINK,
  gsw_average_Res_FarmCPU = gsw_average_Res_FarmCPU,
  gsw_average_Res_MLM = gsw_average_Res_MLM,
  gsw_average_Res_MLMM = gsw_average_Res_MLMM,
  gsw_scaled_sum_Res_BLINK = gsw_scaled_sum_Res_BLINK,
  gsw_scaled_sum_Res_FarmCPU = gsw_scaled_sum_Res_FarmCPU,
  gsw_scaled_sum_Res_MLM = gsw_scaled_sum_Res_MLM,
  gsw_scaled_sum_Res_MLMM = gsw_scaled_sum_Res_MLMM,
  gsw_sum_Res_BLINK = gsw_sum_Res_BLINK,
  gsw_sum_Res_FarmCPU = gsw_sum_Res_FarmCPU,
  gsw_sum_Res_MLM = gsw_sum_Res_MLM,
  gsw_sum_Res_MLMM = gsw_sum_Res_MLMM
)

#Response
df_list <- list(
  response_gsw_diff_first_second_Res_BLINK = response_gsw_diff_first_second_Res_BLINK,
  response_gsw_diff_first_second_Res_FarmCPU = response_gsw_diff_first_second_Res_FarmCPU,
  response_gsw_first_Res_BLINK = response_gsw_first_Res_BLINK,
  response_gsw_first_Res_FarmCPU = response_gsw_first_Res_FarmCPU,
  response_gsw_second_Res_BLINK = response_gsw_second_Res_BLINK,
  response_gsw_second_Res_FarmCPU = response_gsw_second_Res_FarmCPU,
  response_gsw_second_Res_MLMM = response_gsw_second_Res_MLMM,
  response_transp_first_Res_BLINK = response_transp_first_Res_BLINK,
  response_transp_first_Res_FarmCPU = response_transp_first_Res_FarmCPU,
  response_transp_second_Res_BLINK = response_transp_second_Res_BLINK,
  response_transp_second_Res_FarmCPU = response_transp_second_Res_FarmCPU,
  VPDleaf_response_first_Res_BLINK = VPDleaf_response_first_Res_BLINK
)

## Response with bad QQ plots
df_list <- list(
  response_gsw_diff_first_second_Res_BLINK = response_gsw_diff_first_second_Res_BLINK,
  response_gsw_diff_first_second_Res_FarmCPU = response_gsw_diff_first_second_Res_FarmCPU,
  response_gsw_first_Res_BLINK = response_gsw_first_Res_BLINK,
  response_gsw_first_Res_FarmCPU = response_gsw_first_Res_FarmCPU,
  response_gsw_second_Res_BLINK = response_gsw_second_Res_BLINK,
  response_gsw_second_Res_FarmCPU = response_gsw_second_Res_FarmCPU,
  response_gsw_second_Res_MLMM = response_gsw_second_Res_MLMM,
  response_transp_first_Res_BLINK = response_transp_first_Res_BLINK,
  response_transp_first_Res_FarmCPU = response_transp_first_Res_FarmCPU,
  response_transp_second_Res_BLINK = response_transp_second_Res_BLINK,
  response_transp_second_Res_FarmCPU = response_transp_second_Res_FarmCPU,
  VPDleaf_response_first_Res_BLINK = VPDleaf_response_first_Res_BLINK,

  response_gsw_diff_first_second_Res_MLM = response_gsw_diff_first_second_Res_MLM,
  response_gsw_diff_first_second_Res_MLMM = response_gsw_diff_first_second_Res_MLMM,
  response_gsw_first_Res_MLM = response_gsw_first_Res_MLM,
  response_gsw_first_Res_MLMM = response_gsw_first_Res_MLMM,
  response_gsw_second_Res_MLM = response_gsw_second_Res_MLM,
  response_transp_first_Res_MLM = response_transp_first_Res_MLM,
  response_transp_first_Res_MLMM = response_transp_first_Res_MLMM,
  response_transp_second_Res_MLM = response_transp_second_Res_MLM,
  response_transp_second_Res_MLMM = response_transp_second_Res_MLMM,
  VPDleaf_response_first_Res_FarmCPU = VPDleaf_response_first_Res_FarmCPU,
  VPDleaf_response_first_Res_MLM = VPDleaf_response_first_Res_MLM,
  VPDleaf_response_first_Res_MLMM = VPDleaf_response_first_Res_MLMM
)


#Recovery 
df_list <- list(
  recovery_ETR_fifth_Res_BLINK = recovery_ETR_fifth_Res_BLINK,
  recovery_ETR_fifth_Res_MLMM = recovery_ETR_fifth_Res_MLMM,
  recovery_ETR_fifth_Res_FarmCPU = recovery_ETR_fifth_Res_FarmCPU,
  recovery_gsw_diff_fifth_fourth_Res_BLINK = recovery_gsw_diff_fifth_fourth_Res_BLINK,
  recovery_gsw_diff_fifth_fourth_Res_FarmCPU = recovery_gsw_diff_fifth_fourth_Res_FarmCPU,
  recovery_gsw_diff_fifth_fourth_Res_MLMM = recovery_gsw_diff_fifth_fourth_Res_MLMM,
  recovery_gsw_fifth_Res_BLINK = recovery_gsw_fifth_Res_BLINK,
  recovery_gsw_fifth_Res_FarmCPU = recovery_gsw_fifth_Res_FarmCPU,
  recovery_gsw_fifth_Res_MLMM = recovery_gsw_fifth_Res_MLMM,
  recovery_gsw_fourth_Res_BLINK = recovery_gsw_fourth_Res_BLINK,
  recovery_gsw_fourth_Res_FarmCPU = recovery_gsw_fourth_Res_FarmCPU,
  recovery_gsw_fourth_Res_MLMM = recovery_gsw_fourth_Res_MLMM,
  recovery_transp_fourth_Res_BLINK = recovery_transp_fourth_Res_BLINK,
  recovery_transp_fourth_Res_FarmCPU = recovery_transp_fourth_Res_FarmCPU,
  recovery_transp_fourth_Res_MLMM = recovery_transp_fourth_Res_MLMM
)

#Recovery with Bad QQ plots
df_list <- list(
  recovery_ETR_fifth_Res_BLINK = recovery_ETR_fifth_Res_BLINK,
  recovery_ETR_fifth_Res_MLMM = recovery_ETR_fifth_Res_MLMM,
  recovery_ETR_fifth_Res_FarmCPU = recovery_ETR_fifth_Res_FarmCPU,
  recovery_gsw_diff_fifth_fourth_Res_BLINK = recovery_gsw_diff_fifth_fourth_Res_BLINK,
  recovery_gsw_diff_fifth_fourth_Res_FarmCPU = recovery_gsw_diff_fifth_fourth_Res_FarmCPU,
  recovery_gsw_diff_fifth_fourth_Res_MLMM = recovery_gsw_diff_fifth_fourth_Res_MLMM,
  recovery_gsw_fifth_Res_BLINK = recovery_gsw_fifth_Res_BLINK,
  recovery_gsw_fifth_Res_FarmCPU = recovery_gsw_fifth_Res_FarmCPU,
  recovery_gsw_fifth_Res_MLMM = recovery_gsw_fifth_Res_MLMM,
  recovery_gsw_fourth_Res_BLINK = recovery_gsw_fourth_Res_BLINK,
  recovery_gsw_fourth_Res_FarmCPU = recovery_gsw_fourth_Res_FarmCPU,
  recovery_gsw_fourth_Res_MLMM = recovery_gsw_fourth_Res_MLMM,
  recovery_transp_fourth_Res_BLINK = recovery_transp_fourth_Res_BLINK,
  recovery_transp_fourth_Res_FarmCPU = recovery_transp_fourth_Res_FarmCPU,
  recovery_transp_fourth_Res_MLMM = recovery_transp_fourth_Res_MLMM,

  recovery_ETR_fifth_Res_MLM = recovery_ETR_fifth_Res_MLM,
  recovery_gsw_diff_fifth_fourth_Res_MLM = recovery_gsw_diff_fifth_fourth_Res_MLM,
  recovery_gsw_fifth_Res_MLM = recovery_gsw_fifth_Res_MLM,
  recovery_gsw_fourth_Res_MLM = recovery_gsw_fourth_Res_MLM,
  recovery_transp_fourth_Res_MLM = recovery_transp_fourth_Res_MLM
)


#gsw 
df_list <- list(
  gsw_diff_first_fifth_Res_BLINK = gsw_diff_first_fifth_Res_BLINK,
  gsw_diff_first_fifth_Res_FarmCPU = gsw_diff_first_fifth_Res_FarmCPU,
  gsw_diff_first_fifth_Res_MLMM = gsw_diff_first_fifth_Res_MLMM,
  gsw_diff_first_fourth_Res_BLINK = gsw_diff_first_fourth_Res_BLINK,
  gsw_diff_first_fourth_Res_FarmCPU = gsw_diff_first_fourth_Res_FarmCPU,
  gsw_diff_first_fourth_Res_MLMM = gsw_diff_first_fourth_Res_MLMM,
  gsw_diff_first_third_Res_BLINK = gsw_diff_first_third_Res_BLINK,
  gsw_diff_first_third_Res_FarmCPU = gsw_diff_first_third_Res_FarmCPU,
  gsw_diff_first_third_Res_MLMM = gsw_diff_first_third_Res_MLMM,
  gsw_average_Res_BLINK = gsw_average_Res_BLINK,
  gsw_average_Res_FarmCPU = gsw_average_Res_FarmCPU,
  gsw_average_Res_MLM = gsw_average_Res_MLM,
  gsw_average_Res_MLMM = gsw_average_Res_MLMM,
  gsw_scaled_sum_Res_BLINK = gsw_scaled_sum_Res_BLINK,
  gsw_scaled_sum_Res_FarmCPU = gsw_scaled_sum_Res_FarmCPU,
  gsw_scaled_sum_Res_MLM = gsw_scaled_sum_Res_MLM,
  gsw_scaled_sum_Res_MLMM = gsw_scaled_sum_Res_MLMM,
  gsw_sum_Res_BLINK = gsw_sum_Res_BLINK,
  gsw_sum_Res_FarmCPU = gsw_sum_Res_FarmCPU,
  gsw_sum_Res_MLM = gsw_sum_Res_MLM,
  gsw_sum_Res_MLMM = gsw_sum_Res_MLMM
)

#Transp
df_list <- list(
  Transp_average_Res_BLINK = Transp_average_Res_BLINK,
  Transp_average_Res_FarmCPU = Transp_average_Res_FarmCPU,
  Transp_average_Res_MLMM = Transp_average_Res_MLMM,
  Transp_scaled_sum_Res_BLINK = Transp_scaled_sum_Res_BLINK,
  Transp_scaled_sum_Res_FarmCPU = Transp_scaled_sum_Res_FarmCPU,
  Transp_scaled_sum_Res_MLMM = Transp_scaled_sum_Res_MLMM,
  Transp_sum_Res_BLINK = Transp_sum_Res_BLINK,
  Transp_sum_Res_FarmCPU = Transp_sum_Res_FarmCPU,
  Transp_sum_Res_MLMM = Transp_scaled_sum_Res_MLMM,
  Transp_diff_first_fifth_Res_FarmCPU = Transp_diff_first_fifth_Res_FarmCPU,
  Transp_diff_first_fourth_Res_BLINK = Transp_diff_first_fourth_Res_BLINK,
  Transp_diff_first_fourth_Res_FarmCPU = Transp_diff_first_fourth_Res_FarmCPU,
  Transp_diff_first_fourth_Res_MLMM = Transp_diff_first_fourth_Res_MLMM,
  Transp_diff_first_third_Res_BLINK = Transp_diff_first_third_Res_BLINK,
  Transp_diff_first_third_Res_FarmCPU = Transp_diff_first_third_Res_FarmCPU,
  Transp_diff_first_third_Res_MLMM = Transp_diff_first_third_Res_MLMM
)

#Development
df_list <- list(
  No_pods_Res_BLINK = No_pods_Res_BLINK,
  No_pods_Res_FarmCPU = No_pods_Res_FarmCPU,
  Partitioning_weight_Res_BLINK = Partitioning_weight_Res_BLINK,
  Partitioning_weight_Res_FarmCPU = Partitioning_weight_Res_FarmCPU,
  Weight_per_pod_Res_FarmCPU = Weight_per_pod_Res_FarmCPU,
  Weight_pods_Res_FarmCPU = Weight_pods_Res_FarmCPU,
  BBCH_average_Res_BLINK = BBCH_average_Res_BLINK,
  BBCH_average_Res_FarmCPU = BBCH_average_Res_FarmCPU,
  BBCH_diff_fifth_first_Res_BLINK = BBCH_diff_fifth_first_Res_BLINK,
  BBCH_diff_fifth_first_Res_FarmCPU = BBCH_diff_fifth_first_Res_FarmCPU,
  BBCH_diff_fourth_first_Res_BLINK = BBCH_diff_fourth_first_Res_BLINK,
  BBCH_diff_fourth_first_Res_FarmCPU = BBCH_diff_fourth_first_Res_FarmCPU,
  BBCH_diff_third_first_Res_BLINK = BBCH_diff_third_first_Res_BLINK,
  BBCH_sum_Res_BLINK = BBCH_sum_Res_BLINK,
  BBCH_sum_Res_FarmCPU = BBCH_sum_Res_FarmCPU
)

#ETR, VPDleaf and Tleaf
df_list <- list(
  VPDleaf_diff_first_fourth_Res_BLINK = VPDleaf_diff_first_fourth_Res_BLINK,
  VPDleaf_diff_first_fourth_Res_FarmCPU = VPDleaf_diff_first_fourth_Res_FarmCPU,
  VPDleaf_diff_first_fourth_Res_MLMM = VPDleaf_diff_first_fourth_Res_MLMM,
  VPDleaf_diff_first_third_Res_MLMM = VPDleaf_diff_first_third_Res_MLMM,
  Tleaf_diff_first_fourth_Res_FarmCPU = Tleaf_first_fourth_Res_FarmCPU,
  Tleaf_diff_first_fourth_Res_MLMM = Tleaf_first_fourth_Res_MLMM,
  Tleaf_diff_first_third_Res_MLMM = Tleaf_first_third_Res_MLMM,
  ETR_diff_first_fifth_Res_BLINK = ETR_diff_first_fifth_Res_BLINK,
  ETR_diff_first_fifth_Res_FarmCPU = ETR_diff_first_fifth_Res_FarmCPU,
  ETR_diff_first_fifth_Res_MLM = ETR_diff_first_fifth_Res_MLM,
  ETR_diff_first_fifth_Res_MLMM = ETR_diff_first_fifth_Res_MLMM
)

##Thinning dfs for plotting

thinned_list <- list()

for (i in seq_along(df_list)) {
  df <- df_list[[i]]
  
  # Split into high and low significance
  high_sig <- df %>% filter(-log10(P.value) > 3)
  low_sig  <- df %>% filter(-log10(P.value) <= 3)
  
  # Sample a fraction (e.g., 10%) from the low significance group
  low_sig_sampled <- sample_frac(low_sig, 0.1)
  
  # Combine and store
  thinned_list[[i]] <- bind_rows(high_sig, low_sig_sampled)
}

names(thinned_list) <- names(df_list)

processed_list <- list()

for (df_name in names(thinned_list)) {
  x <- thinned_list[[df_name]]
  
  # Extract last part of the name
  model_name <- sub(".*_(\\w+)$", "\\1", df_name)
  
  #Extract trait name
  trait_name <- sub("_Res_.*", "", df_name)
  
  # Compute cumulative position
  don <- x %>%
    group_by(Chr) %>%
    summarise(chr_len = max(Pos), .groups = 'drop') %>%
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    left_join(x, ., by = "Chr") %>%
    arrange(Chr, Pos) %>%
    mutate(BPcum = Pos + tot)
  
  # Add model column
  don$model <- model_name
  don$trait <- trait_name
  
  # Store in list
  processed_list[[df_name]] <- don
}

don <- processed_list[[1]]

axisdf = don %>%
 group_by(Chr) %>%
 summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

gwasResults <- do.call(rbind, processed_list)
rownames(gwasResults) <- NULL
gwasResults$Trait_model <- paste(gwasResults$trait, gwasResults$model, sep = "_")

gwasResults$Sizes <- -log10(gwasResults$P.value)
gwasResults$size <- scales::rescale(gwasResults$Sizes, to = c(1, 7))
gwasResults$Observed <- -log10(gwasResults$P.value)


sig_all <- -log10(0.05 / nrow(don))

gwasResults$Trait_model <- gsub( '_', ' ', gwasResults$Trait_model)
gwasResults$Trait_model <- gsub( 'fifth', 'W5', gwasResults$Trait_model)
gwasResults$Trait_model <- gsub( 'fourth', 'W4', gwasResults$Trait_model)
gwasResults$Trait_model <- gsub( 'third', 'W3', gwasResults$Trait_model)
gwasResults$Trait_model <- gsub( 'second', 'W2', gwasResults$Trait_model)
gwasResults$Trait_model <- gsub( 'first', 'W1', gwasResults$Trait_model)
gwasResults$trait <- gsub( '_', ' ', gwasResults$trait)
gwasResults$trait <- gsub( 'fifth', 'W5', gwasResults$trait)
gwasResults$trait <- gsub( 'fourth', 'W4', gwasResults$trait)
gwasResults$trait <- gsub( 'third', 'W3', gwasResults$trait)
gwasResults$trait <- gsub( 'second', 'W2', gwasResults$trait)
gwasResults$trait <- gsub( 'first', 'W1', gwasResults$trait)

gwasResults$Trait_model <- gsub( 'gsw', 'Gsw', gwasResults$Trait_model)
gwasResults$trait <- gsub( 'gsw', 'Gsw', gwasResults$trait)

gwasResults$Trait_model <- gsub( 'No pods', 'Pod number', gwasResults$Trait_model)
gwasResults$Trait_model <- gsub( 'Weight pods', 'Pod weight', gwasResults$Trait_model)
gwasResults$trait <- gsub( 'No pods', 'Pod number', gwasResults$trait)
gwasResults$trait <- gsub( 'Weight pods', 'Pod weight', gwasResults$trait)

gwasResults$Trait_model <- gsub( 'transp', 'E', gwasResults$Trait_model)
gwasResults$Trait_model <- gsub( 'Transp', 'E', gwasResults$Trait_model)
gwasResults$trait <- gsub( 'transp', 'E', gwasResults$trait)
gwasResults$trait <- gsub( 'Transp', 'E', gwasResults$trait)

gwasResults$Trait_model <- gsub( 'response', 'Response', gwasResults$Trait_model)
gwasResults$trait <- gsub( 'response', 'Response', gwasResults$trait)

gwasResults$Trait_model <- gsub( 'recovery', 'Recovery', gwasResults$Trait_model)
gwasResults$trait <- gsub( 'recovery', 'Recovery', gwasResults$trait)

#gwasResults$Trait_model <- sub("^ETR\\s+", "", gwasResults$Trait_model)

#### Vertical lines --------------------------------------------
##Add vertical lines for matching
new <- gwasResults[gwasResults$Chr == "Chr11", ] 
new2 <- new[new$Pos > 8043419,]

#gsw
#Chr06:24894074, Chr11:6147062, Chr11:8043420
ver_lines <- c(268405931, 466535353, 468431711)

#gsw diff
#Chr01:6780993, Chr02:35721884, Chr04:510279, Chr04:9525571, Chr04:41296047, 
#Chr06:3327906, Chr08:9402343, Chr08:45430325, Chr11:2271453, Chr11:8043420
ver_lines <- c(6780993, 87155707, 155053236, 164068528, 195839004, 
               246839763, 324191177, 360219159, 462659744, 468431711)

#gsw recovery 
#Chr01:6780993, Chr02:5443201, Chr02:35721884, Chr03:4619356, Chr03:11060166, 
#Chr04:41272563, Chr04:41296047, Chr05:17841150, Chr05:38027262, Chr07:28259272, 
#Chr08:16621250, Chr08:61914250, Chr09:29402247, Chr10:12228269, 
#Chr10:35434949, Chr11:6147062, Chr11:8043420
ver_lines <- c(6780993, 56877024, 87155707, 105723781, 112164591, 195815520, 195839004,
               220432450, 240618562, 303007186, 331410084, 376703084, 407238671, 
               428314708, 451521388, 466535353, 468431711
               )

#gsw response
#Chr04:510279, Chr05:39947225, Chr06:3327906, Chr06:27468233,
#Chr08:9402343, Chr08:45430325, Chr11:24147900
ver_lines <- c(155053236, 242538525, 246839763, 270980090, 324191177,
               360219159, 484536191
                )

#Transp only
#Chr06:3327906, Chr06:24894074, Chr09:29402247
ver_lines <- c(246839763, 268405931, 407238671)

#Transp diff
#Chr02:27176203, Chr02:45713765, Chr05:39947225, Chr06:3327906, Chr08:9402343
#Chr08:45430325
ver_lines <- c(78610026, 97147588, 242538525, 246839763, 324191177, 360219159)

#Transp recovery and response
#Chr02:5443201, Chr02:27178537, Chr03:4619356, Chr03:11060166, Chr04:41272563
#Chr05:38027262, Chr06:27468233, Chr08:9402343, Chr08:16621250, Chr08:61914250
#Chr09:29402247, Chr10:12228269, Chr10:35434949, Chr11:24147900
ver_lines <- c(56877024, 78612360, 105723781, 112164591, 195815520, 240618562,
               270980090, 324191177, 331410084, 376703084, 407238671, 428314708,
               451521388, 484536191
               )

#BBCH vertical lines
#Chr04:9523515, Chr07:650547, Chr07:28251449, Chr08:343107, Chr08:51734523,
#Chr09:7557718, Chr11:1007337
ver_lines <- c(164066472, 275398461, 302999363, 315131941, 366523357,
               385394142, 461395628)

#Harvesting vertical lines
#Chr02:45713765, Chr07:650547, Chr08:343107, Chr08:51734523, Chr09:7557718, 
#Chr11:1007337
ver_lines <- c(97147588, 275398461, 315131941, 366523357, 385394142, 461395628)

#Tleaf vertical lines
#Chr11:12365041, Chr01:25238058
ver_lines <- c(472753332, 25238058)

#VPDleaf vertical lines
#Chr11:12365041, Chr01:25238058
ver_lines <- c(472753332, 25238058)

#ETR vertial lines
#Chr11: 16821066
ver_lines <- c(477209357)

#Response vertical lines
ver_lines <- c(78612360, 155053236, 242538525, 246839763, 270980090, 324191177,
               360219159, 484536191)

#Recovery vertical lines
ver_lines <- c(6780993, 56877024, 87155707, 105723781,  112164591, 
               195815520, 195839004, 220432450, 240618562, 303007186, 
               331410084, 376703084, 407238671, 428314708, 451521388, 
               466535353, 468431711, 477209357)

#gsw vertical lines
ver_lines <- c(6780993, 87155707, 155053236, 164068528, 195839004, 
               246839763,268405931, 324191177, 360219159, 462659744, 
               466535353, 468431711)

#Transp vertical lines
ver_lines <- c(78610026, 97147588, 242538525, 246839763, 268405931, 
               324191177, 360219159, 407238671)

#BBCH, harvesting
ver_lines <- c(164066472, 275398461, 302999363, 315131941, 366523357,
               385394142, 461395628, 97147588)

#ETR, VPDleaf, Tleaf
ver_lines <- c(472753332, 25238058, 477209357)

#Check
#-log10(1.977080e-36)
x <- -log10(min(gwasResults$P.value))
max <- ceiling(x / 10) * 10


### Plotting all --------------------------------------------
shape_values <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)
"CHANGE P"
p5 <- ggplot(gwasResults, aes(x=BPcum, y=-log10(P.value), shape = Trait_model)) +
  scale_shape_manual(values = shape_values) + 
  # Show all points
  geom_point( aes(color=as.factor(Chr)), alpha=1, size=gwasResults$size, stroke = 1.5) +
  #scale_color_manual(values = rep(c("#ff6666", "#ffbd55", "#9de24f", "#87cefa", "#cf9dff"), 11 )) +
  scale_color_manual(values = rep(c("#87cefa", "darkgrey"), 11 )) +
  
  
  #Bonferonni cutoff
  geom_hline(
    yintercept = sig_all, color = "forestgreen", linewidth = 1.4) +
  
  ##Add vertical lines for matching with Andean or Whole panel  
  geom_vline(xintercept = ver_lines, color = "black", linetype = "dotted", linewidth = 1.4) +
  
  
  # custom X axis:
  scale_x_continuous(label = axisdf$Chr, breaks= axisdf$center, expand = c(0.01, 0) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max)  ) +     # remove space between plot area and x axis
  
  guides(colour = "none") +
  
  # Custom the theme:
  theme_classic() +
  theme( 
    axis.title.x=element_blank(),
    legend.position="bottom",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size=18),
    axis.title.y = element_text(size=18),
    legend.text = element_text(size=16),
    legend.title = element_text(size=18)
  ) + guides(shape = guide_legend(title = "Trait and model", 
                                  override.aes = list(size = 5), nrow = 3)) 

#p, p1, p2, p3, p4, p5
#gsw, transp, development, ETR/VPD/Tleaf, response, recovery

#pdf(file ="det_A_gwas_F_manhattan.pdf", width = 10, height = 4)
#plot(p)
#dev.off()

#jpeg("det_A_gwas_MF_manhattan.jpg", width = 1500, height = 600)
"CHANGE P"
png("Response_leaf_manhattan_gaps_higher_y_smaller_newnames.png",  width = 1800, height = 600)
plot(p4)
dev.off()

#p 

### ggarrange
library(ggpubr)
library(png)
library(grid)
library(ggplot2)

# Function to wrap a PNG image in a ggplot of fixed dimensions
# image_to_plot <- function(path) {
#   img <- rasterGrob(readPNG(path), interpolate = TRUE)
#   
#   ggplot() +
#     annotation_custom(img, xmin = -1, xmax = 1, ymin = -1, ymax = 1) +
#     coord_fixed() +
#     theme_void()
# }
# 
# p1 <- image_to_plot("gsw/gsw_all_manhattan_gaps_higher_y.png")
# p2 <- image_to_plot("gsw/gsw_diff_manhattan_gaps_higher_y.png")
# p3 <- image_to_plot("gsw/gsw_recovery_manhattan_gaps_higher_y.png")
# p4 <- image_to_plot("gsw/gsw_response_manhattan_gaps_higher_y.png")

# p1 <- rasterGrob(readPNG("gsw/gsw_all_manhattan_gaps_higher_y.png"), interpolate = TRUE)
# p2 <- rasterGrob(readPNG("gsw/gsw_diff_manhattan_gaps_higher_y.png"), interpolate = TRUE)
# p3 <- rasterGrob(readPNG("gsw/gsw_recovery_manhattan_gaps_higher_y.png"), interpolate = TRUE)
# p4 <- rasterGrob(readPNG("gsw/gsw_response_manhattan_gaps_higher_y.png"), interpolate = TRUE)
# 
# p1 <- p1 + theme(legend.position = "bottom")
# p2 <- p2 + theme(legend.position = "bottom")
# p3 <- p3 + theme(legend.position = "bottom")
# p4 <- p4 + theme(legend.position = "bottom")

png("manhattan__all_newnames_2.png",  width = 1600, height = 2400)
print(ggarrange(
  p, p1, p2, p3, p4, p5,
  ncol = 1, 
  nrow = 6, 
  labels = "AUTO", 
  font.label = list(size = 24)
))
dev.off()


##### QQ plots ----------------------------------
#install.packages("ggbreak")
library(ggbreak)
#library(ggforce)

df_subset <- gwasResults %>%
  arrange(P.value) %>%              # sort by the column of interest
  mutate(row_id = row_number()) %>%     # add row number
  filter(row_id %% 17 == 1) %>%         # select every 10th row (starting at 1st)
  select(-row_id)  

write.csv(df_subset, "QQ_plots/thinned_df_subset_recovery.csv")

#df_subset <- read.csv("QQ_plots/thinned_df_subset_gsw.csv")
df_subset <- read.csv("QQ_plots/thinned_df_subset_transp.csv")
#df_subset <- read.csv("QQ_plots/thinned_df_subset_harvesting.csv")
#df_subset <- read.csv("QQ_plots/thinned_df_subset_VET.csv") #VET #response #recovery
#df_subset <- read.csv("QQ_plots/thinned_df_subset_response.csv") #VET #response #recovery
#df_subset <- read.csv("QQ_plots/thinned_df_subset_recovery.csv") #VET #response #recovery


x <- -log10(min(df_subset$P.value))
max <- ceiling(x / 10) * 10

colnames(df_subset)[13] <- c("Trait")

safe_pal_8 <- c("#000000", "#888888", "#5EC962", "#117733", "#332288", "#6699CC", "#CC5500","#AA4499")


"CHANGE QQ"
qq2 <- ggplot(df_subset) +
  stat_qq(size = 1, aes(sample=Observed, colour=Trait)) +
  geom_abline(slope=1, intercept=0) +
  #scale_colour_brewer(palette = "Dark2") +
  scale_fill_manual(values=safe_pal_8) +
  xlab("Expected -log10(p)") +
  ylab("Observed -log10(p)") +
  
  scale_x_continuous(expand = c(0, 0), limits = c(0,5) ) +
  #scale_y_continuous(limits = c(0, 16)) +   #16 For VET
  #scale_y_continuous(limits = c(0, 35)) +   #For response
  scale_y_continuous(limits = c(0, max)) +
  
  # Custom the theme:
  theme_gray() +
  theme( 
    axis.text = element_text(size=10),
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size=12),
    legend.text = element_text(size=10),
    legend.title = element_text(size=12), 
    axis.text.y.right  = element_blank(),
    axis.ticks.y.right = element_blank()) +
  facet_wrap(~model, ncol = 4, axes = "all") +
  theme(panel.spacing.x = unit(1, "lines")) +
  #scale_y_break(c(16, 35), ticklabels = c(35, 40), space = 0.1) #GSW
  scale_y_break(c(13, 18), ticklabels = c(18, 20), space = 0.1) #transp
  #scale_y_break(c(16, 45), ticklabels = c(45, 50), space = 0.1) #harvesting
  #scale_y_break(c(15, 30), ticklabels = c(30, 35), space = 0.2) #response
  #scale_y_break(c(25, 45), ticklabels = c(45, 50), space = 0.2) + #recovery
  #scale_y_break(c(50, 95), ticklabels = c(95, 100), space = 0.2) #recovery
  
#GSW, Transp, Harvesting, VET, Response, Recovery
#qq1, qq2, qq3, qq4, qq5, qq6

"CHANGE QQ"
# pdf(file ="QQ_plots/all_QQ_1.pdf", width = 9, height = 4)
# qq1
# dev.off()

ggsave("QQ_plots/all_QQ_6.pdf", plot = {print(qq6)}, width = 8.5, height = 3.5)

qq1_m <- qq1 + theme(plot.margin = margin(t = 0.2, r= 2, b = 0.2, l = 0.2, "cm"))
qq2_m <- qq2 + theme(plot.margin = margin(t = 0.2, r= 2.5, b = 0.2, l = 0.2, "cm"))
qq3_m <- qq3 + theme(plot.margin = margin(t = 0.2, r= 1.6, b = 0.2, l = 0.2, "cm"))
qq4_m <- qq4 + theme(plot.margin = margin(t = 0.2, r= 1.3, b = 0.2, l = 0.2, "cm"))


library(ggpubr)

#pdf(file ="QQ_plots/All_QQ_new.pdf", width = 9, height = 13)
gg_all <- ggarrange(
  print(qq1), print(qq2), print(qq3), print(qq4),nrow = 4, ncol = 1,  
  labels = "AUTO"
)
gg_all
#dev.off()

gg_all <- ggarrange(
  print(qq5), print(qq6),nrow = 2, ncol = 1, 
  labels = c("E", "F")
)

ggsave("QQ_plots/all_QQ_break.pdf", plot = {print(gg_all)}, width = 9, height = 14)

