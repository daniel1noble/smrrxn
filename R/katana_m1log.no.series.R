setwd("~/gitrepo/smrrxn/")

#clear envir
rm(list = ls())

#load library
require(MCMCglmm)
require(parallel)

#set seed
set.seed(1)

#read in data

data <- read.csv("data/data_final/mrrxn_final_v2.csv")
data$id <- as.factor(data$id)
data$series <- as.factor(data$series)

varibs.need <- c("obs", "samp_period", "id" ,"batch", "series", "incb_num", "incb_temp_id", "defecate", "incb_temp", "z.incb_temp", "z.log.temp", "incb_temp_K", "inverseK_incb_temp", "body_temp", "z.body_temp", "z.log.body_temp", "body_temp_K", "inverseK_body_temp" , "z.prior_temp1", "z.log.prior_temp1", "prior_temp1_K", "inverseK_prior_temp1", "prior_temp2_K", "inverseK_prior_temp2", "z.prior_temp2", "z.log.prior_temp2", "orig_lizmass", "lizmass_nocombout", "log.mass", "z.log.mass","orig_co2_pmin", "co2pm_nocombout", "log.co2pmin", "z.log.co2pmin")

incl.vars <- names(data) %in% varibs.need
data <- data[incl.vars]

predictors <- c("log.co2pmin", "inverseK_incb_temp", "z.log.mass", "inverseK_prior_temp2", "id", "series")

dat <- data[complete.cases(data[,predictors]),]
str(dat)

#priors
expanded.prior <- list(R = list(V = 1, nu = 0.002),
                       G = list(G1 = list(V = diag(2), nu = 0.002, alpha.V = diag(1000,2,2), alpha.mu = rep(0,2))))


#m1.log.noprior - removes prior temp and series
m1.log.noprior.noseries <- mclapply(1:3, function(i) {
  MCMCglmm(log.co2pmin ~ inverseK_incb_temp + z.log.mass + samp_period,
          random = ~us(1+inverseK_incb_temp):id,
          family = "gaussian",
          prior = expanded.prior,
           nitt = 7510000,
           burnin = 10000,
           thin = 5000,
           data = dat, 
           verbose = T)
}, mc.cores = 3)

saveRDS(m1.log.noprior.noseries, "R/m1.log.noprior.noseries")



