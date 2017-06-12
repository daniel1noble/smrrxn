setwd("~/Dropbox/smrrxn/")

rm(list = ls())

library(lme4)

rawdata <- read.csv("data/data_final/ldeli_smr.csv")
str(rawdata)

#standardising variables

data <- rawdata[,1:17]

co2_mlpmin ~ incb_temp + lizmass + (1+incb_temp|id) + (1+incb_temp|series_id) 




