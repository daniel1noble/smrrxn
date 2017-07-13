setwd("~/Dropbox/smrrxn/")

#clear envir
rm(list = ls())

#load library
library(MCMCglmm)

#read in data

data <- read.csv("data/data_final/mrrxn_final.csv")
data$id <- as.factor(data$id)
data$series <- as.factor(data$series)

varibs.need <- c("obs", "samp_period", "id" ,"batch", "series", "incb_num", "incb_temp_id", "defecate", "incb_temp", "z.incb_temp", "z.log.temp", "incb_temp_K", "inverseK_incb_temp", "body_temp", "z.body_temp", "z.log.body_temp", "body_temp_K", "inverseK_body_temp" , "z.prior_temp1", "z.log.prior_temp1", "prior_temp1_K", "inverseK_prior_temp1", "prior_temp2_K", "inverseK_prior_temp2", "z.prior_temp2", "z.log.prior_temp2", "orig_lizmass", "lizmass_nocombout", "log.mass", "z.log.mass","orig_co2_pmin", "co2pm_nocombout", "log.co2pmin", "z.log.co2pmin")

incl.vars <- names(data) %in% varibs.need
data <- data[incl.vars]

predictors <- c("z.log.co2pmin", "inverseK_incb_temp", "z.log.mass", "inverseK_prior_temp2", "id", "series")

dat <- data[complete.cases(data[,predictors]),]
str(dat)

#priors
expanded.prior <- list(R = list(V = 1, nu = 0.002),
                       G = list(G1 = list(V = diag(2), nu = 0.002, alpha.V = diag(1000,2,2), alpha.mu = rep(0,2)),
                                G2 = list(V = diag(2), nu = 0.002, alpha.V = diag(1000,2,2), alpha.mu = rep(0,2))))

IW.prior <- list(R = list(V = 1, nu = 0.002),
                    G = list(G1 = list(V = diag(2), nu = 0.002),
                             G2 = list(V = diag(2), nu = 0.002)))

#final model - 1 chain for now
model.1 <- MCMCglmm(z.log.co2pmin ~ inverseK_incb_temp + z.log.mass + inverseK_prior_temp2,
                    random = ~us(1+inverseK_incb_temp):id + us(1+inverseK_incb_temp):series,
                    family = "gaussian",
                    prior = expanded.prior,
                    nitt = 5010000,
                    burnin = 10000,
                    thin = 5000,
                    data = dat, 
                    verbose = T)


saveRDS(model.1, "output/model.1")

plot(model.1)
summary(model.1)


