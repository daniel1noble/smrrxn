setwd("~/gitrepo/smrrxn/")

#load library
library(MCMCglmm)
library(parallel)

data <- read.csv("data/data_final/mr_final_recentered.csv")
data$id <- as.factor(data$id)
data$series <- as.factor(data$series)

#priors
expanded.prior <- list(R = list(V = 1, nu = 0.002),
                       G = list(G1 = list(V = diag(2), nu = 0.002, alpha.V = diag(1000,2,2), alpha.mu = rep(0,2)),
                                G2 = list(V = diag(2), nu = 0.002, alpha.V = diag(1000,2,2), alpha.mu = rep(0,2))))


m1log_t22 <- FALSE
m1log_t24 <- FALSE
m1log_t26 <- FALSE
m1log_t28 <- TRUE
m1log_t30 <- FALSE
m1log_t32 <- FALSE


#22 centered
if(m1log_t22){
  m1log_t22 <- mclapply(1:3, function(i) {
    MCMCglmm(log.co2pmin ~ inverseK_incb_temp_22cen + z.log.mass,
             random = ~us(1+inverseK_incb_temp_22cen):id + us(1+inverseK_incb_temp_22cen):series,
             family = "gaussian",
             prior = expanded.prior,
             nitt = 7510000,
             burnin = 10000,
             thin = 5000,
             data = data, 
             verbose = T)
  }, mc.cores = 3)
  saveRDS(m1log_t22, "R/m1log_t22")
}

#24 centered
if(m1log_t24){
  m1log_t24 <- mclapply(1:3, function(i) {
    MCMCglmm(log.co2pmin ~ inverseK_incb_temp_24cen + z.log.mass,
             random = ~us(1+inverseK_incb_temp_24cen):id + us(1+inverseK_incb_temp_24cen):series,
             family = "gaussian",
             prior = expanded.prior,
             nitt = 7510000,
             burnin = 10000,
             thin = 5000,
             data = data, 
             verbose = T)
  }, mc.cores = 3)
  saveRDS(m1log_t24, "R/m1log_t24")
}

#26 centered
if(m1log_t26){
  m1log_t26 <- mclapply(1:3, function(i) {
    MCMCglmm(log.co2pmin ~ inverseK_incb_temp_26cen + z.log.mass,
             random = ~us(1+inverseK_incb_temp_26cen):id + us(1+inverseK_incb_temp_26cen):series,
             family = "gaussian",
             prior = expanded.prior,
             nitt = 7510000,
             burnin = 10000,
             thin = 5000,
             data = data, 
             verbose = T)
  }, mc.cores = 3)
  saveRDS(m1log_t26, "R/m1log_t26")
}


#28 centered
if(m1log_t28){
  m1log_t28 <- mclapply(1:3, function(i) {
    MCMCglmm(log.co2pmin ~ inverseK_incb_temp_28cen + z.log.mass,
             random = ~us(1+inverseK_incb_temp_28cen):id + us(1+inverseK_incb_temp_28cen):series,
             family = "gaussian",
             prior = expanded.prior,
             nitt = 7510000,
             burnin = 10000,
             thin = 5000,
             data = data, 
             verbose = T)
  }, mc.cores = 3)
  saveRDS(m1log_t28, "R/m1log_t28")
}

#30 centered
if(m1log_t30){
  m1log_t30 <- mclapply(1:3, function(i) {
    MCMCglmm(log.co2pmin ~ inverseK_incb_temp_30cen + z.log.mass,
             random = ~us(1+inverseK_incb_temp_30cen):id + us(1+inverseK_incb_temp_30cen):series,
             family = "gaussian",
             prior = expanded.prior,
             nitt = 7510000,
             burnin = 10000,
             thin = 5000,
             data = data, 
             verbose = T)
  }, mc.cores = 3)
  saveRDS(m1log_t30, "R/m1log_t30")
}

#32 centered
if(m1log_t32){
  m1log_t32 <- mclapply(1:3, function(i) {
    MCMCglmm(log.co2pmin ~ inverseK_incb_temp_32cen + z.log.mass,
             random = ~us(1+inverseK_incb_temp_32cen):id + us(1+inverseK_incb_temp_32cen):series,
             family = "gaussian",
             prior = expanded.prior,
             nitt = 7510000,
             burnin = 10000,
             thin = 5000,
             data = data, 
             verbose = T)
  }, mc.cores = 3)
  saveRDS(m1log_t32, "R/m1log_t32")
}