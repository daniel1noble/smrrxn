setwd("~/gitrepo/smrrxn/")

#clear envir
rm(list = ls())

#load library
library(MCMCglmm)

#read in data

data <- read.csv("data/data_final/mr_final_analysis.csv")

#priors

expanded.prior <- list(R = list(V = 1, nu = 0.002),
                       G = list(G1 = list(V = diag(2), nu = 0.002, alpha.V = rep(1000,2,2), alpha.mu = rep(0,2)),
                                G2 = list(V = diag(2), nu = 0.002, alpha.V = rep(1000,2,2), alpha.mu = rep(0,2))))

#final model - 1 chain for now

model.1 <- MCMCglmm(z.log.co2pmin ~ z.log.temp + z.log.mass + z.prior_temp_4,
                    random = ~us(1+z.log.temp):id + ~us(1+trial):series,
                    family = "gaussian",
                    prior = my.prior,
                    nitt = 5010000,
                    burnin = 10000,
                    thin = 5000,
                    data = data)

saveRDS(model.1, file="~/output/model.1")

