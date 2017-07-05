setwd("~/Dropbox/smrrxn/")

rm(list = ls())

library(ggplot2)
library(lme4)





write.csv(data, "data/data_final/mr_data_analysis.csv", row.names = F)


data <- read.csv("data/data_final/mr_data_analysis_example.csv")

#standardising variables

data$z.incb_temp <- scale(data$incb_temp)

data$z.body_temp <- scale(data$body_temp)

data$log.mass <- log(data$lizmass)
data$z.log.mass <- scale(data$log.mass)

data$log.co2pmin <- log(data$co2_pmin)
data$z.log.co2pmin <- scale(data$log.co2pmin)

#Creating prior temp variables 
data$z.prior_temp_3 <- scale(data$prior_temp_3)
data$z.prior_temp_4 <- scale(data$prior_temp_4)

data$z.log.prior_temp_3 <- scale(log(data$prior_temp_3))
data$z.log.prior_temp_4 <- scale(log(data$prior_temp_4))

#Log temp
data$z.log.temp <- scale(log(data$incb_temp))

write.csv(data, row.names = F, "data/data_final/mr_final_analysis.csv")





