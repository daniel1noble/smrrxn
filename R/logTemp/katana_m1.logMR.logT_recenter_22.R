setwd("~/gitrepo/smrrxn/")
#setwd("~/Dropbox/smrrxn/")

#load library
library(MCMCglmm)
library(parallel)

data <- read.csv("data/data_final/mr_final_log.T_recentered.csv")
data$id <- as.factor(data$id)
data$series <- as.factor(data$series)

predictors <- c("log.temp_22cen", "z.log.mass", "id", "series", "log.prior_temp1")

data <- data[complete.cases(data[,predictors]),]

#priors
expanded.prior <- list(R = list(V = 1, nu = 0.002),
                       G = list(G1 = list(V = diag(2), nu = 0.002, alpha.V = diag(1000,2,2), alpha.mu = rep(0,2)),
                                G2 = list(V = diag(2), nu = 0.002, alpha.V = diag(1000,2,2), alpha.mu = rep(0,2))))

#22 centered
m1.logMR.logT.t22 <- mclapply(1:3, function(i) {
  MCMCglmm(log.co2pmin ~ log.temp_22cen + z.log.mass + log.prior_temp1,
           random = ~us(1+log.temp_22cen):id + us(1+log.temp_22cen):series,
           family = "gaussian",
           prior = expanded.prior,
           nitt = 7510000,
           burnin = 10000,
           thin = 5000,
           data = data, 
           verbose = T)
}, mc.cores = 3)
  
saveRDS(m1.logMR.logT.t22, "R/m1.logMR.logT.t22")
  


  

