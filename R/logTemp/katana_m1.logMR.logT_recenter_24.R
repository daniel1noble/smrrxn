setwd("~/gitrepo/smrrxn/")

#load library
library(MCMCglmm)
library(parallel)

data <- read.csv("data/data_final/mr_final_log.T_recentered.csv")
data$id <- as.factor(data$id)
data$series <- as.factor(data$series)

predictors <- c("log.temp_24cen", "z.log.mass", "id", "series")

data <- data[complete.cases(data[,predictors]),]

#priors
expanded.prior <- list(R = list(V = 1, nu = 0.002),
                       G = list(G1 = list(V = diag(2), nu = 0.002, alpha.V = diag(1000,2,2), alpha.mu = rep(0,2)),
                                G2 = list(V = diag(2), nu = 0.002, alpha.V = diag(1000,2,2), alpha.mu = rep(0,2))))

#24 centered
m1.logMR.logT.npT.t24 <- mclapply(1:3, function(i) {
    MCMCglmm(log.co2pmin ~ log.temp_24cen + z.log.mass,
             random = ~us(1+log.temp_24cen):id + us(1+log.temp_24cen):series,
             family = "gaussian",
             prior = expanded.prior,
             nitt = 7510000,
             burnin = 10000,
             thin = 5000,
             data = data, 
             verbose = T)
  }, mc.cores = 3)

  saveRDS(m1.logMR.logT.npT.t24, "R/m1.logMR.logT.npT.t24")


