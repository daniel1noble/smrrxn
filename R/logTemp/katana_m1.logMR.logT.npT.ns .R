setwd("~/gitrepo/smrrxn/")

#clear envir
rm(list = ls())

#load library
require(MCMCglmm)
require(parallel)

#set seed
set.seed(1)

#read in data

data <- read.csv("data/data_final/mrrxn_logT.csv")
data$id <- as.factor(data$id)
data$series <- as.factor(data$series)

varibs.need <- c("samp_period", "id" , "series", "log.temp", "log.prior_temp1",  "log.prior_temp2", "z.log.mass","log.co2pmin", "z.log.co2pmin")

incl.vars <- names(data) %in% varibs.need
data <- data[incl.vars]

predictors <- c("log.co2pmin", "log.temp", "z.log.mass", "id", "series", "samp_period")

dat <- data[complete.cases(data[,predictors]),]
str(dat)

#priors
expanded.prior <- list(R = list(V = 1, nu = 0.002),
                       G = list(G1 = list(V = diag(2), nu = 0.002, alpha.V = diag(1000,2,2), alpha.mu = rep(0,2))))


#m1.logMR.logT.npT.ns - removes prior temp and series
m1.logMR.logT.npT.ns <- mclapply(1:3, function(i) {
  MCMCglmm(log.co2pmin ~ log.temp + z.log.mass + samp_period,
          random = ~us(1+log.temp):id,
          family = "gaussian",
          prior = expanded.prior,
           nitt = 7510000,
           burnin = 10000,
           thin = 5000,
           data = dat, 
           verbose = T)
}, mc.cores = 3)

saveRDS(m1.logMR.logT.npT.ns, "R/m1.logMR.logT.npT.ns")



