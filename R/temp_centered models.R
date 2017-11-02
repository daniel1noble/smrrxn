setwd("~/Dropbox/smrrxn/")

#clear envir
rm(list = ls())

#load library
library(MCMCglmm)
library(parallel)
library(dplyr)
library(ggplot2)

source("R/functions/smr_functions.R")

####################################
#### INTERCEPT at 22 degrees C #####
####################################

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

##Rslope
t22_R.slope <- m1_t22.VCV[,"inverseK_incb_temp_22cen:inverseK_incb_temp_22cen.id"] / ( m1_t22.VCV[,"inverseK_incb_temp_22cen:inverseK_incb_temp_22cen.id"] + m1_t22.VCV[,"inverseK_incb_temp_22cen:inverseK_incb_temp_22cen.series"] ) 

posterior.mode(t22_R.slope) #Rslope
HPDinterval(as.mcmc(t22_R.slope)) #Rslope CIs

##Covariance/Correlation of intercept and slope ID level

m1_t22.id.cor <- m1_t22.VCV[,"inverseK_incb_temp_22cen:(Intercept).id"] / ( sqrt(m1_t22.VCV[,"(Intercept):(Intercept).id"]) * sqrt(m1_t22.VCV[,"inverseK_incb_temp_22cen:inverseK_incb_temp_22cen.id"]) )
posterior.mode(m1_t22_cor)
HPDinterval(as.mcmc(m1_t22_cor))

##Covariance/Correlation of intercept and slope series level

m1_t22.series.cor <- m1_t22.VCV[,"inverseK_incb_temp_22cen:(Intercept).series"] / ( sqrt(m1_t22.VCV[,"(Intercept):(Intercept).series"]) * sqrt(m1_t22.VCV[,"inverseK_incb_temp_22cen:inverseK_incb_temp_22cen.series"]) )
posterior.mode(m1_t22.series.cor)
HPDinterval(as.mcmc(m1_t22.series.cor))

####################################
#### INTERCEPT at 24 degrees C #####
####################################

m1_t24 <- readRDS("output/rds/m1_t24")

#Diagnostics
m1_t24.S <- lapply(m1_t24, function(m) m$Sol)
m1_t24.Sol <- do.call(rbind, m1_t24.S)
m1_t24.S1 <- do.call(mcmc.list, m1_t24.S)

m1_t24.V <- lapply(m1_t24, function(m) m$VCV)
m1_t24.VCV <- do.call(rbind, m1_t24.V)
m1_t24.V1 <- do.call(mcmc.list, m1_t24.V)

plot(m1_t24.S1)
gelman.diag(m1_t24.S1)
summary(m1_t24.Sol)
posterior.mode(m1_t24.Sol) ; HPDinterval(as.mcmc(m1_t24.Sol))

plot(m1_t24.V1)
gelman.diag(m1_t24.V1, multivariate = F)
summary(m1_t24.VCV)
posterior.mode(m1_t24.VCV) ; HPDinterval(as.mcmc(m1_t24.VCV))

##Rint
t24_R.int <- m1_t24.VCV[,"(Intercept):(Intercept).id"] / ( m1_t24.VCV[,"(Intercept):(Intercept).id"] + m1_t24.VCV[,"(Intercept):(Intercept).series"] ) 

posterior.mode(t24_R.int) #R.int
HPDinterval(as.mcmc(t24_R.int)) #R.int CIs

##Rslope
t24_R.slope <- m1_t24.VCV[,"inverseK_incb_temp_24cen:inverseK_incb_temp_24cen.id"] / ( m1_t24.VCV[,"inverseK_incb_temp_24cen:inverseK_incb_temp_24cen.id"] + m1_t24.VCV[,"inverseK_incb_temp_24cen:inverseK_incb_temp_24cen.series"] ) 

posterior.mode(t24_R.slope) #Rslope
HPDinterval(as.mcmc(t24_R.slope)) #Rslope CIs

##Covariance/Correlation of intercept and slope ID level

m1_t24.id.cor <- m1_t24.VCV[,"inverseK_incb_temp_24cen:(Intercept).id"] / ( sqrt(m1_t24.VCV[,"(Intercept):(Intercept).id"]) * sqrt(m1_t24.VCV[,"inverseK_incb_temp_24cen:inverseK_incb_temp_24cen.id"]) )
posterior.mode(m1_t24_cor)
HPDinterval(as.mcmc(m1_t24_cor))

##Covariance/Correlation of intercept and slope series level

m1_t24.series.cor <- m1_t24.VCV[,"inverseK_incb_temp_24cen:(Intercept).series"] / ( sqrt(m1_t24.VCV[,"(Intercept):(Intercept).series"]) * sqrt(m1_t24.VCV[,"inverseK_incb_temp_24cen:inverseK_incb_temp_24cen.series"]) )
posterior.mode(m1_t24.series.cor)
HPDinterval(as.mcmc(m1_t24.series.cor))

####################################
#### INTERCEPT at 26 degrees C #####
####################################

m1_t26 <- readRDS("output/rds/m1_t26")

#Diagnostics
m1_t26.S <- lapply(m1_t26, function(m) m$Sol)
m1_t26.Sol <- do.call(rbind, m1_t26.S)
m1_t26.S1 <- do.call(mcmc.list, m1_t26.S)

m1_t26.V <- lapply(m1_t26, function(m) m$VCV)
m1_t26.VCV <- do.call(rbind, m1_t26.V)
m1_t26.V1 <- do.call(mcmc.list, m1_t26.V)

plot(m1_t26.S1)
gelman.diag(m1_t26.S1)
summary(m1_t26.Sol)
posterior.mode(m1_t26.Sol) ; HPDinterval(as.mcmc(m1_t26.Sol))

plot(m1_t26.V1)
gelman.diag(m1_t26.V1, multivariate = F)
summary(m1_t26.VCV)
posterior.mode(m1_t26.VCV) ; HPDinterval(as.mcmc(m1_t26.VCV))

##Rint
t26_R.int <- m1_t26.VCV[,"(Intercept):(Intercept).id"] / ( m1_t26.VCV[,"(Intercept):(Intercept).id"] + m1_t26.VCV[,"(Intercept):(Intercept).series"] ) 

posterior.mode(t26_R.int) #R.int
HPDinterval(as.mcmc(t26_R.int)) #R.int CIs

##Rslope
t26_R.slope <- m1_t26.VCV[,"inverseK_incb_temp_26cen:inverseK_incb_temp_26cen.id"] / ( m1_t26.VCV[,"inverseK_incb_temp_26cen:inverseK_incb_temp_26cen.id"] + m1_t26.VCV[,"inverseK_incb_temp_26cen:inverseK_incb_temp_26cen.series"] ) 

posterior.mode(t26_R.slope) #Rslope
HPDinterval(as.mcmc(t26_R.slope)) #Rslope CIs

##Covariance/Correlation of intercept and slope ID level

m1_t26.id.cor <- m1_t26.VCV[,"inverseK_incb_temp_26cen:(Intercept).id"] / ( sqrt(m1_t26.VCV[,"(Intercept):(Intercept).id"]) * sqrt(m1_t26.VCV[,"inverseK_incb_temp_26cen:inverseK_incb_temp_26cen.id"]) )
posterior.mode(m1_t26_cor)
HPDinterval(as.mcmc(m1_t26_cor))

##Covariance/Correlation of intercept and slope series level

m1_t26.series.cor <- m1_t26.VCV[,"inverseK_incb_temp_26cen:(Intercept).series"] / ( sqrt(m1_t26.VCV[,"(Intercept):(Intercept).series"]) * sqrt(m1_t26.VCV[,"inverseK_incb_temp_26cen:inverseK_incb_temp_26cen.series"]) )
posterior.mode(m1_t26.series.cor)
HPDinterval(as.mcmc(m1_t26.series.cor))

####################################
#### INTERCEPT at 28 degrees C #####
####################################

m1_t28 <- readRDS("output/rds/m1_t28")

#Diagnostics
m1_t28.S <- lapply(m1_t28, function(m) m$Sol)
m1_t28.Sol <- do.call(rbind, m1_t28.S)
m1_t28.S1 <- do.call(mcmc.list, m1_t28.S)

m1_t28.V <- lapply(m1_t28, function(m) m$VCV)
m1_t28.VCV <- do.call(rbind, m1_t28.V)
m1_t28.V1 <- do.call(mcmc.list, m1_t28.V)

plot(m1_t28.S1)
gelman.diag(m1_t28.S1)
summary(m1_t28.Sol)
posterior.mode(m1_t28.Sol) ; HPDinterval(as.mcmc(m1_t28.Sol))

plot(m1_t28.V1)
gelman.diag(m1_t28.V1, multivariate = F)
summary(m1_t28.VCV)
posterior.mode(m1_t28.VCV) ; HPDinterval(as.mcmc(m1_t28.VCV))

##Rint
t28_R.int <- m1_t28.VCV[,"(Intercept):(Intercept).id"] / ( m1_t28.VCV[,"(Intercept):(Intercept).id"] + m1_t28.VCV[,"(Intercept):(Intercept).series"] ) 

posterior.mode(t28_R.int) #R.int
HPDinterval(as.mcmc(t28_R.int)) #R.int CIs

##Rslope
t28_R.slope <- m1_t28.VCV[,"inverseK_incb_temp_28cen:inverseK_incb_temp_28cen.id"] / ( m1_t28.VCV[,"inverseK_incb_temp_28cen:inverseK_incb_temp_28cen.id"] + m1_t28.VCV[,"inverseK_incb_temp_28cen:inverseK_incb_temp_28cen.series"] ) 

posterior.mode(t28_R.slope) #Rslope
HPDinterval(as.mcmc(t28_R.slope)) #Rslope CIs

##Covariance/Correlation of intercept and slope ID level

m1_t28.id.cor <- m1_t28.VCV[,"inverseK_incb_temp_28cen:(Intercept).id"] / ( sqrt(m1_t28.VCV[,"(Intercept):(Intercept).id"]) * sqrt(m1_t28.VCV[,"inverseK_incb_temp_28cen:inverseK_incb_temp_28cen.id"]) )
posterior.mode(m1_t28_cor)
HPDinterval(as.mcmc(m1_t28_cor))

##Covariance/Correlation of intercept and slope series level

m1_t28.series.cor <- m1_t28.VCV[,"inverseK_incb_temp_28cen:(Intercept).series"] / ( sqrt(m1_t28.VCV[,"(Intercept):(Intercept).series"]) * sqrt(m1_t28.VCV[,"inverseK_incb_temp_28cen:inverseK_incb_temp_28cen.series"]) )
posterior.mode(m1_t28.series.cor)
HPDinterval(as.mcmc(m1_t28.series.cor))

####################################
#### INTERCEPT at 30 degrees C #####
####################################

m1_t30 <- readRDS("output/rds/m1_t30")

#Diagnostics
m1_t30.S <- lapply(m1_t30, function(m) m$Sol)
m1_t30.Sol <- do.call(rbind, m1_t30.S)
m1_t30.S1 <- do.call(mcmc.list, m1_t30.S)

m1_t30.V <- lapply(m1_t30, function(m) m$VCV)
m1_t30.VCV <- do.call(rbind, m1_t30.V)
m1_t30.V1 <- do.call(mcmc.list, m1_t30.V)

plot(m1_t30.S1)
gelman.diag(m1_t30.S1)
summary(m1_t30.Sol)
posterior.mode(m1_t30.Sol) ; HPDinterval(as.mcmc(m1_t30.Sol))

plot(m1_t30.V1)
gelman.diag(m1_t30.V1, multivariate = F)
summary(m1_t30.VCV)
posterior.mode(m1_t30.VCV) ; HPDinterval(as.mcmc(m1_t30.VCV))

##Rint
t30_R.int <- m1_t30.VCV[,"(Intercept):(Intercept).id"] / ( m1_t30.VCV[,"(Intercept):(Intercept).id"] + m1_t30.VCV[,"(Intercept):(Intercept).series"] ) 

posterior.mode(t30_R.int) #R.int
HPDinterval(as.mcmc(t30_R.int)) #R.int CIs

##Rslope
t30_R.slope <- m1_t30.VCV[,"inverseK_incb_temp_30cen:inverseK_incb_temp_30cen.id"] / ( m1_t30.VCV[,"inverseK_incb_temp_30cen:inverseK_incb_temp_30cen.id"] + m1_t30.VCV[,"inverseK_incb_temp_30cen:inverseK_incb_temp_30cen.series"] ) 

posterior.mode(t30_R.slope) #Rslope
HPDinterval(as.mcmc(t30_R.slope)) #Rslope CIs

##Covariance/Correlation of intercept and slope ID level

m1_t30.id.cor <- m1_t30.VCV[,"inverseK_incb_temp_30cen:(Intercept).id"] / ( sqrt(m1_t30.VCV[,"(Intercept):(Intercept).id"]) * sqrt(m1_t30.VCV[,"inverseK_incb_temp_30cen:inverseK_incb_temp_30cen.id"]) )
posterior.mode(m1_t30_cor)
HPDinterval(as.mcmc(m1_t30_cor))

##Covariance/Correlation of intercept and slope series level

m1_t30.series.cor <- m1_t30.VCV[,"inverseK_incb_temp_30cen:(Intercept).series"] / ( sqrt(m1_t30.VCV[,"(Intercept):(Intercept).series"]) * sqrt(m1_t30.VCV[,"inverseK_incb_temp_30cen:inverseK_incb_temp_30cen.series"]) )
posterior.mode(m1_t30.series.cor)
HPDinterval(as.mcmc(m1_t30.series.cor))

####################################
#### INTERCEPT at 32 degrees C #####
####################################

m1_t32 <- readRDS("output/rds/m1_t32")

#Diagnostics
m1_t32.S <- lapply(m1_t32, function(m) m$Sol)
m1_t32.Sol <- do.call(rbind, m1_t32.S)
m1_t32.S1 <- do.call(mcmc.list, m1_t32.S)

m1_t32.V <- lapply(m1_t32, function(m) m$VCV)
m1_t32.VCV <- do.call(rbind, m1_t32.V)
m1_t32.V1 <- do.call(mcmc.list, m1_t32.V)

plot(m1_t32.S1)
gelman.diag(m1_t32.S1)
summary(m1_t32.Sol)
posterior.mode(m1_t32.Sol) ; HPDinterval(as.mcmc(m1_t32.Sol))

plot(m1_t32.V1)
gelman.diag(m1_t32.V1, multivariate = F)
summary(m1_t32.VCV)
posterior.mode(m1_t32.VCV) ; HPDinterval(as.mcmc(m1_t32.VCV))

##Rint
t32_R.int <- m1_t32.VCV[,"(Intercept):(Intercept).id"] / ( m1_t32.VCV[,"(Intercept):(Intercept).id"] + m1_t32.VCV[,"(Intercept):(Intercept).series"] ) 

posterior.mode(t32_R.int) #R.int
HPDinterval(as.mcmc(t32_R.int)) #R.int CIs

##Rslope
t32_R.slope <- m1_t32.VCV[,"inverseK_incb_temp_32cen:inverseK_incb_temp_32cen.id"] / ( m1_t32.VCV[,"inverseK_incb_temp_32cen:inverseK_incb_temp_32cen.id"] + m1_t32.VCV[,"inverseK_incb_temp_32cen:inverseK_incb_temp_32cen.series"] ) 

posterior.mode(t32_R.slope) #Rslope
HPDinterval(as.mcmc(t32_R.slope)) #Rslope CIs

##Covariance/Correlation of intercept and slope ID level

m1_t32.id.cor <- m1_t32.VCV[,"inverseK_incb_temp_32cen:(Intercept).id"] / ( sqrt(m1_t32.VCV[,"(Intercept):(Intercept).id"]) * sqrt(m1_t32.VCV[,"inverseK_incb_temp_32cen:inverseK_incb_temp_32cen.id"]) )
posterior.mode(m1_t32_cor)
HPDinterval(as.mcmc(m1_t32_cor))

##Covariance/Correlation of intercept and slope series level

m1_t32.series.cor <- m1_t32.VCV[,"inverseK_incb_temp_32cen:(Intercept).series"] / ( sqrt(m1_t32.VCV[,"(Intercept):(Intercept).series"]) * sqrt(m1_t32.VCV[,"inverseK_incb_temp_32cen:inverseK_incb_temp_32cen.series"]) )
posterior.mode(m1_t32.series.cor)
HPDinterval(as.mcmc(m1_t32.series.cor))



