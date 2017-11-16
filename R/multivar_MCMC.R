setwd("~/gitrepo/smrrxn/")

#clear envir
rm(list = ls())

#load library
library(MCMCglmm)
library(parallel)

#data
data <- read.csv("data/data_final/multiresp_log.csv")
data$id <- as.factor(data$id)
data$series <- as.factor(data$series)

#varibs.need <- c("obs", "samp_period", "id" ,"batch", "series", "incb_num", "incb_temp_id", "defecate", "incb_temp", "z.incb_temp", "z.log.temp", "incb_temp_K", "inverseK_incb_temp", "body_temp", "z.body_temp", "z.log.body_temp", "body_temp_K", "inverseK_body_temp" , "z.prior_temp1", "z.log.prior_temp1", "prior_temp1_K", "inverseK_prior_temp1", "prior_temp2_K", "inverseK_prior_temp2", "z.prior_temp2", "z.log.prior_temp2", "orig_lizmass", "lizmass_nocombout","orig_lizmass", "log.mass", "z.log.mass","orig_co2_pmin", "co2pm_nocombout", "log.co2pmin", "z.log.co2pmin")

#incl.vars <- names(data) %in% varibs.need
#data <- data[incl.vars]

#data3 <- data %>% 
  #select(id, samp_period, series, incb_temp, log.co2pmin, orig_lizmass) %>%
  #group_by(series) %>%
  #mutate(avg_mass = mean(orig_lizmass)) %>%
  #select(id, samp_period, series, incb_temp, log.co2pmin, avg_mass) %>%
  #spread(incb_temp, log.co2pmin) 
  
#colnames(data3)[5:10] <- paste0("t_", colnames(data3)[5:10])

#data3$z.log.mass <- scale(log(data3$avg_mass))

#data3 <- as.data.frame(data3)
#data3 <- data3[!is.na(data3$z.log.mass),]
#write.csv(data3, row.names = F, "data/data_final/multiresp_log.csv")

#priors
multi.prior <- list(R = list(V = diag(6), nu = 0.01), G = list(G1 = list(V = diag(6), nu = 0.01)))

#m4
 m4_sp <- mclapply(1:3, function(i) {
  MCMCglmm(cbind(t_22, t_24, t_26, t_28, t_30, t_32) ~  z.log.mass + samp_period,
             random= ~us(trait):id,
             rcov = ~us(trait):units,
             family = c(rep("gaussian", 6)),
             prior = multi.prior,
             nitt = 7510000,
             burnin = 10000,
             thin = 5000,
             data = data, 
             verbose = T)
  }, mc.cores = 3)
  
 saveRDS(m4_sp, "R/m4_sp_noz")






