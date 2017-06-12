setwd("~/Dropbox/smrrxn/")

rm(list = ls())

library(ggplot2)
library(lme4)

rawdata <- read.csv("data/data_final/ldelipop_smr.csv")
str(rawdata)

rawdata <- rawdata[order(rawdata$samp_period),]

data <- rawdata[,1:18]
#names(data)
[1] "date"         "samp_period"  "series_id"    "batch"        "incb_num"     "incb_temp_id"
[7] "incb_temp"    "id"           "tc"           "body_temp"    "lizmass"      "ch_vol"      
[13] "total_air"    "t_diff"       "frac_co2"     "total_co2"    "co2_mlpmin"   "defacate"   

summary(data$co2_mlpmin)




#standardising variables








#Analysis
co2_mlpmin ~ incb_temp + lizmass + (1+incb_temp|id) + (1+incb_temp|series_id) 




