setwd("~/Dropbox/smrrxn/")

#clear envir
rm(list = ls())

#load library
library(MCMCglmm)
library(parallel)
library(tidyr)
library(dplyr)

#data
data <- read.csv("data/data_final/mrrxn_final_v2.csv")
data$id <- as.factor(data$id)
data$series <- as.factor(data$series)

head(data)
names(data)

varibs.need <- c("obs", "samp_period", "id" ,"batch", "series", "incb_num", "incb_temp_id", "defecate", "incb_temp", "z.incb_temp", "z.log.temp", "incb_temp_K", "inverseK_incb_temp", "body_temp", "z.body_temp", "z.log.body_temp", "body_temp_K", "inverseK_body_temp" , "z.prior_temp1", "z.log.prior_temp1", "prior_temp1_K", "inverseK_prior_temp1", "prior_temp2_K", "inverseK_prior_temp2", "z.prior_temp2", "z.log.prior_temp2", "orig_lizmass", "lizmass_nocombout","orig_lizmass", "log.mass", "z.log.mass","orig_co2_pmin", "co2pm_nocombout", "log.co2pmin", "z.log.co2pmin")

incl.vars <- names(data) %in% varibs.need
data <- data[incl.vars]

head(data)
names(data)

data3 <- data %>% 
  select(id, samp_period, series, incb_temp, z.log.co2pmin, orig_lizmass) %>%
  group_by(series) %>%
  mutate(avg_mass = mean(orig_lizmass)) %>%
  select(id, samp_period, series, incb_temp, z.log.co2pmin, avg_mass) %>%
  spread(incb_temp, z.log.co2pmin)

colnames(data3)[5:10] <- paste0("t_", colnames(data3)[5:10])
View(data3)

write.csv(data3, row.names = F, "data/data_final/multirandom.csv")






