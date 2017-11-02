setwd("~/Dropbox/smrrxn/")

#clear envir
rm(list = ls())

#load library
library(MCMCglmm)
library(parallel)
library(dplyr)
library(ggplot2)

source("R/functions/smr_functions.R")

#### INTERCEPT at 22 degrees C #####
m1_t22 <- readRDS("output/rds/m1_t22")

#Diagnostics
m1_t22.S <- lapply(m1_t22, function(m) m$Sol)
m1_t22.Sol <- do.call(rbind, m1_t22.S)
m1_t22.S1 <- do.call(mcmc.list, m1_t22.S)

m1_t22.V <- lapply(m1_t22, function(m) m$VCV)
m1_t22.VCV <- do.call(rbind, m1_t22.V)
m1_t22.V1 <- do.call(mcmc.list, m1_t22.V)

plot(m1_t22.S1)
gelman.diag(m1_t22.S1)
summary(m1_t22.Sol)
posterior.mode(m1_t22.Sol) ; HPDinterval(as.mcmc(m1_t22.Sol))

plot(m1_t22.V1)
gelman.diag(m1_t22.V1, multivariate = F)
summary(m1_t22.VCV)
posterior.mode(m1_t22.VCV) ; HPDinterval(as.mcmc(m1_t22.VCV))

##Rint
t22_R.int <- m1_t22.VCV[,"(Intercept):(Intercept).id"] / ( m1_t22.VCV[,"(Intercept):(Intercept).id"] + m1_t22.VCV[,"(Intercept):(Intercept).series"] ) 

posterior.mode(t22_R.int) #R.int
HPDinterval(as.mcmc(t22_R.int)) #R.int CIs

##Covariance/Correlation of intercept and slope ID level

m1_t22.id.cor <- m1_t22.VCV[,"inverseK_incb_temp_22cen:(Intercept).id"] / ( sqrt(m1_t22.VCV[,"(Intercept):(Intercept).id"]) * sqrt(m1_t22.VCV[,"inverseK_incb_temp_22cen:inverseK_incb_temp_22cen.id"]) )
posterior.mode(m1_t22_cor)
HPDinterval(as.mcmc(m1_t22_cor))

##Covariance/Correlation of intercept and slope series level

m1_t22.series.cor <- m1_t22.VCV[,"inverseK_incb_temp_22cen:(Intercept).series"] / ( sqrt(m1_t22.VCV[,"(Intercept):(Intercept).series"]) * sqrt(m1_t22.VCV[,"inverseK_incb_temp_22cen:inverseK_incb_temp_22cen.series"]) )
posterior.mode(m1_t22.series.cor)
HPDinterval(as.mcmc(m1_t22.series.cor))

#### INTERCEPT at 24 degrees C #####