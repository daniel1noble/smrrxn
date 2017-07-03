setwd("~/gitrepo/smrrxn/")

#clear envir
rm(list = ls())

#load library
library(MCMCglmm)

#read in data

data <- read.csv("data/data_final/mr_final_analysis.csv")

#priors

my.prior <- list(V =1, nu = 0.002, alpha.V = 1000, alpha.mu = 0)

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

