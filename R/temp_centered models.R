setwd("~/Dropbox/smrrxn/")

#clear envir
rm(list = ls())

#load library
library(MCMCglmm)
library(parallel)
library(dplyr)
library(ggplot2)

source("R/functions/smr_functions.R")

#####Setting up tables for data######
TableRint <- data.frame(matrix(nrow = 6 , ncol = 9))
rownames(TableRint) <- c(22, 24, 26, 28, 30, 32)
colnames(TableRint) <- c("R_int", "R_lower", "R_upper",
                         "ID_var", "ID_lower", "ID_upper",
                         "Series_var", "Series_lower", "Series_upper")
TableRint$Temperature <- rownames(TableRint)

TableRslope <- data.frame(matrix(nrow = 6 , ncol = 9))
rownames(TableRslope) <- c(22, 24, 26, 28, 30, 32)
colnames(TableRslope) <- c("R_slope", "R_lower", "R_upper",
                           "ID_var", "ID_lower", "ID_upper",
                           "Series_var", "Series_lower", "Series_upper")
TableRslope$Temperature <- rownames(TableRslope)

TableRshort <- data.frame(matrix(nrow = 6 , ncol = 3))
rownames(TableRshort) <- c(22, 24, 26, 28, 30, 32)
colnames(TableRshort) <- c("Repeatability", "lower", "upper")

TableRlong <- data.frame(matrix(nrow = 6 , ncol = 3))
rownames(TableRlong) <- c(22, 24, 26, 28, 30, 32)
colnames(TableRlong) <- c("Repeatability", "lower", "upper")

TableIntSlopeCovCor <- data.frame(matrix(nrow = 6 , ncol = 6))
rownames(TableIntSlopeCovCor) <- c(22, 24, 26, 28, 30, 32)
colnames(TableIntSlopeCovCor) <- c("Covariance", "Cov_lower", "Cov_upper", "Correlation", "Cor_lower", "Cor_upper")

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

TableRint[1,1] <- posterior.mode(t22_R.int) #R.int
TableRint[1,2:3] <-HPDinterval(as.mcmc(t22_R.int)) #R.int CIs

#Rint Between ID variance
TableRint[1,4] <- posterior.mode(m1_t22.VCV[,"(Intercept):(Intercept).id"]) 
TableRint[1,5:6] <- HPDinterval(as.mcmc(m1_t22.VCV[,"(Intercept):(Intercept).id"]))

#Rint Series ID (Within ID - amoung series) variance
TableRint[1,7] <- posterior.mode(m1_t22.VCV[,"(Intercept):(Intercept).series"]) 
TableRint[1,8:9] <- HPDinterval(as.mcmc(m1_t22.VCV[,"(Intercept):(Intercept).series"]))

##Rslope
t22_R.slope <- m1_t22.VCV[,"inverseK_incb_temp_22cen:inverseK_incb_temp_22cen.id"] / ( m1_t22.VCV[,"inverseK_incb_temp_22cen:inverseK_incb_temp_22cen.id"] + m1_t22.VCV[,"inverseK_incb_temp_22cen:inverseK_incb_temp_22cen.series"] ) 

TableRslope[1,1] <- posterior.mode(t22_R.slope) #Rslope
TableRslope[1,2:3] <- HPDinterval(as.mcmc(t22_R.slope)) #Rslope CIs

#Rslope between ID variance
TableRslope[1,4] <- posterior.mode(m1_t22.VCV[,"inverseK_incb_temp_22cen:inverseK_incb_temp_22cen.id"])
TableRslope[1,5:6] <- HPDinterval(as.mcmc(m1_t22.VCV[,"inverseK_incb_temp_22cen:inverseK_incb_temp_22cen.id"]))

#Rslope within ID variance amoung series variance
TableRslope[1,7] <- posterior.mode(m1_t22.VCV[,"inverseK_incb_temp_22cen:inverseK_incb_temp_22cen.series"])
TableRslope[1,8:9] <- HPDinterval(as.mcmc(m1_t22.VCV[,"inverseK_incb_temp_22cen:inverseK_incb_temp_22cen.series"]))

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

TableRint[2,1] <- posterior.mode(t24_R.int) #R.int
TableRint[2,2:3] <-HPDinterval(as.mcmc(t24_R.int)) #R.int CIs

#Between ID variance
TableRint[2,4] <- posterior.mode(m1_t24.VCV[,"(Intercept):(Intercept).id"]) 
TableRint[2,5:6] <- HPDinterval(as.mcmc(m1_t24.VCV[,"(Intercept):(Intercept).id"]))

#Series ID (Within ID - amoung series) variance
TableRint[2,7] <- posterior.mode(m1_t24.VCV[,"(Intercept):(Intercept).series"]) 
TableRint[2,8:9] <- HPDinterval(as.mcmc(m1_t24.VCV[,"(Intercept):(Intercept).series"]))

##Rslope
t24_R.slope <- m1_t24.VCV[,"inverseK_incb_temp_24cen:inverseK_incb_temp_24cen.id"] / ( m1_t24.VCV[,"inverseK_incb_temp_24cen:inverseK_incb_temp_24cen.id"] + m1_t24.VCV[,"inverseK_incb_temp_24cen:inverseK_incb_temp_24cen.series"] ) 

TableRslope[2,1] <- posterior.mode(t24_R.slope) #Rslope
TableRslope[2,2:3] <- HPDinterval(as.mcmc(t24_R.slope)) #Rslope CIs

#Rslope between ID variance
TableRslope[2,4] <- posterior.mode(m1_t24.VCV[,"inverseK_incb_temp_24cen:inverseK_incb_temp_24cen.id"])
TableRslope[2,5:6] <- HPDinterval(as.mcmc(m1_t24.VCV[,"inverseK_incb_temp_24cen:inverseK_incb_temp_24cen.id"]))

#Rslope within ID variance amoung series variance
TableRslope[2,7] <- posterior.mode(m1_t24.VCV[,"inverseK_incb_temp_24cen:inverseK_incb_temp_24cen.series"])
TableRslope[2,8:9] <- HPDinterval(as.mcmc(m1_t24.VCV[,"inverseK_incb_temp_24cen:inverseK_incb_temp_24cen.series"]))

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

TableRint[3,1] <- posterior.mode(t26_R.int)  #R.int
TableRint[3,2:3] <-HPDinterval(as.mcmc(t26_R.int)) #R.int CIs

#Between ID variance
TableRint[3,4] <- posterior.mode(m1_t26.VCV[,"(Intercept):(Intercept).id"]) 
TableRint[3,5:6] <- HPDinterval(as.mcmc(m1_t26.VCV[,"(Intercept):(Intercept).id"]))

#Series ID (Within ID - amoung series) variance
TableRint[3,7] <- posterior.mode(m1_t26.VCV[,"(Intercept):(Intercept).series"]) 
TableRint[3,8:9] <- HPDinterval(as.mcmc(m1_t26.VCV[,"(Intercept):(Intercept).series"]))

##Rslope
t26_R.slope <- m1_t26.VCV[,"inverseK_incb_temp_26cen:inverseK_incb_temp_26cen.id"] / ( m1_t26.VCV[,"inverseK_incb_temp_26cen:inverseK_incb_temp_26cen.id"] + m1_t26.VCV[,"inverseK_incb_temp_26cen:inverseK_incb_temp_26cen.series"] ) 

TableRslope[3,1] <- posterior.mode(t26_R.slope) #Rslope
TableRslope[3,2:3] <- HPDinterval(as.mcmc(t26_R.slope)) #Rslope CIs

#Rslope between ID variance
TableRslope[3,4] <- posterior.mode(m1_t26.VCV[,"inverseK_incb_temp_26cen:inverseK_incb_temp_26cen.id"])
TableRslope[3,5:6] <- HPDinterval(as.mcmc(m1_t26.VCV[,"inverseK_incb_temp_26cen:inverseK_incb_temp_26cen.id"]))

#Rslope within ID variance amoung series variance
TableRslope[3,7] <- posterior.mode(m1_t26.VCV[,"inverseK_incb_temp_26cen:inverseK_incb_temp_26cen.series"])
TableRslope[3,8:9] <- HPDinterval(as.mcmc(m1_t26.VCV[,"inverseK_incb_temp_26cen:inverseK_incb_temp_26cen.series"]))

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

TableRint[4,1] <- posterior.mode(t28_R.int) #R.int
TableRint[4,2:3] <- HPDinterval(as.mcmc(t28_R.int)) #R.int CIs

#Between ID variance
TableRint[4,4] <- posterior.mode(m1_t28.VCV[,"(Intercept):(Intercept).id"]) 
TableRint[4,5:6] <- HPDinterval(as.mcmc(m1_t28.VCV[,"(Intercept):(Intercept).id"]))

#Series ID (Within ID - amoung series) variance
TableRint[4,7] <- posterior.mode(m1_t28.VCV[,"(Intercept):(Intercept).series"]) 
TableRint[4,8:9] <- HPDinterval(as.mcmc(m1_t28.VCV[,"(Intercept):(Intercept).series"]))

##Rslope
t28_R.slope <- m1_t28.VCV[,"inverseK_incb_temp_28cen:inverseK_incb_temp_28cen.id"] / ( m1_t28.VCV[,"inverseK_incb_temp_28cen:inverseK_incb_temp_28cen.id"] + m1_t28.VCV[,"inverseK_incb_temp_28cen:inverseK_incb_temp_28cen.series"] ) 

TableRslope[4,1] <- posterior.mode(t28_R.slope) #Rslope
TableRslope[4,2:3] <- HPDinterval(as.mcmc(t28_R.slope)) #Rslope CIs

#Rslope between ID variance
TableRslope[4,4] <- posterior.mode(m1_t28.VCV[,"inverseK_incb_temp_28cen:inverseK_incb_temp_28cen.id"])
TableRslope[4,5:6] <- HPDinterval(as.mcmc(m1_t28.VCV[,"inverseK_incb_temp_28cen:inverseK_incb_temp_28cen.id"]))

#Rslope within ID variance amoung series variance
TableRslope[4,7] <- posterior.mode(m1_t28.VCV[,"inverseK_incb_temp_28cen:inverseK_incb_temp_28cen.series"])
TableRslope[4,8:9] <- HPDinterval(as.mcmc(m1_t28.VCV[,"inverseK_incb_temp_28cen:inverseK_incb_temp_28cen.series"]))

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

TableRint[5,1] <- posterior.mode(t30_R.int) #R.int
TableRint[5,2:3] <- HPDinterval(as.mcmc(t30_R.int)) #R.int CIs

#Between ID variance
TableRint[5,4] <- posterior.mode(m1_t30.VCV[,"(Intercept):(Intercept).id"]) 
TableRint[5,5:6] <- HPDinterval(as.mcmc(m1_t30.VCV[,"(Intercept):(Intercept).id"]))

#Series ID (Within ID - amoung series) variance
TableRint[5,7] <- posterior.mode(m1_t30.VCV[,"(Intercept):(Intercept).series"]) 
TableRint[5,8:9] <- HPDinterval(as.mcmc(m1_t30.VCV[,"(Intercept):(Intercept).series"]))

##Rslope
t30_R.slope <- m1_t30.VCV[,"inverseK_incb_temp_30cen:inverseK_incb_temp_30cen.id"] / ( m1_t30.VCV[,"inverseK_incb_temp_30cen:inverseK_incb_temp_30cen.id"] + m1_t30.VCV[,"inverseK_incb_temp_30cen:inverseK_incb_temp_30cen.series"] ) 

TableRslope[5,1] <- posterior.mode(t30_R.slope) #Rslope
TableRslope[5,2:3] <- HPDinterval(as.mcmc(t30_R.slope)) #Rslope CIs

#Rslope between ID variance
TableRslope[5,4] <- posterior.mode(m1_t30.VCV[,"inverseK_incb_temp_30cen:inverseK_incb_temp_30cen.id"])
TableRslope[5,5:6] <- HPDinterval(as.mcmc(m1_t30.VCV[,"inverseK_incb_temp_30cen:inverseK_incb_temp_30cen.id"]))

#Rslope within ID variance amoung series variance
TableRslope[5,7] <- posterior.mode(m1_t30.VCV[,"inverseK_incb_temp_30cen:inverseK_incb_temp_30cen.series"])
TableRslope[5,8:9] <- HPDinterval(as.mcmc(m1_t30.VCV[,"inverseK_incb_temp_30cen:inverseK_incb_temp_30cen.series"]))

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

TableRint[6,1] <- posterior.mode(t32_R.int) #R.int
TableRint[6,2:3] <- HPDinterval(as.mcmc(t32_R.int)) #R.int CIs

#Between ID variance
TableRint[6,4] <- posterior.mode(m1_t32.VCV[,"(Intercept):(Intercept).id"]) 
TableRint[6,5:6] <- HPDinterval(as.mcmc(m1_t32.VCV[,"(Intercept):(Intercept).id"]))

#Series ID (Within ID - amoung series) variance
TableRint[6,7] <- posterior.mode(m1_t32.VCV[,"(Intercept):(Intercept).series"]) 
TableRint[6,8:9] <- HPDinterval(as.mcmc(m1_t32.VCV[,"(Intercept):(Intercept).series"]))

round(TableRint[1:8],2)

write.csv(round(TableRint[1:8],2), row.names = F, "output/table/FV_Rint_varcomp.csv")

##Rslope
t32_R.slope <- m1_t32.VCV[,"inverseK_incb_temp_32cen:inverseK_incb_temp_32cen.id"] / ( m1_t32.VCV[,"inverseK_incb_temp_32cen:inverseK_incb_temp_32cen.id"] + m1_t32.VCV[,"inverseK_incb_temp_32cen:inverseK_incb_temp_32cen.series"] ) 

TableRslope[6,1] <- posterior.mode(t32_R.slope) #Rslope
TableRslope[6,2:3] <- HPDinterval(as.mcmc(t32_R.slope)) #Rslope CIs

#Rslope between ID variance
TableRslope[6,4] <- posterior.mode(m1_t32.VCV[,"inverseK_incb_temp_32cen:inverseK_incb_temp_32cen.id"])
TableRslope[6,5:6] <- HPDinterval(as.mcmc(m1_t32.VCV[,"inverseK_incb_temp_32cen:inverseK_incb_temp_32cen.id"]))

#Rslope within ID variance amoung series variance
TableRslope[6,7] <- posterior.mode(m1_t32.VCV[,"inverseK_incb_temp_32cen:inverseK_incb_temp_32cen.series"])
TableRslope[6,8:9] <- HPDinterval(as.mcmc(m1_t32.VCV[,"inverseK_incb_temp_32cen:inverseK_incb_temp_32cen.series"]))

write.csv(round(TableRslope[1:8],2), row.names = F, "output/table/FV_Rslope_varcomp.csv")

##Covariance/Correlation of intercept and slope ID level

m1_t32.id.cor <- m1_t32.VCV[,"inverseK_incb_temp_32cen:(Intercept).id"] / ( sqrt(m1_t32.VCV[,"(Intercept):(Intercept).id"]) * sqrt(m1_t32.VCV[,"inverseK_incb_temp_32cen:inverseK_incb_temp_32cen.id"]) )
posterior.mode(m1_t32_cor)
HPDinterval(as.mcmc(m1_t32_cor))

##Covariance/Correlation of intercept and slope series level

m1_t32.series.cor <- m1_t32.VCV[,"inverseK_incb_temp_32cen:(Intercept).series"] / ( sqrt(m1_t32.VCV[,"(Intercept):(Intercept).series"]) * sqrt(m1_t32.VCV[,"inverseK_incb_temp_32cen:inverseK_incb_temp_32cen.series"]) )
posterior.mode(m1_t32.series.cor)
HPDinterval(as.mcmc(m1_t32.series.cor))

posterior.mode(m1_t32.series.cor) #R.int
HPDinterval(as.mcmc(m1_t32.series.cor)) #R.int CIs

#Plotting FV forest plots - Repeatability in Intercepts
Rint_FV1 <- ggplot(data = TableRint, aes(x = Temperature, y = Repeatability)) +
  geom_point() + 
  geom_errorbar(aes(ymin = R_lower, ymax = R_upper), width = 0) +
  scale_y_continuous(limits = c(0,1)) + 
  #scale_x_continuous(breaks = x) + 
  labs(y = "Repeatability", x = "Temperature") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

#Plotting FV forest plots - Between ID variance in intercept
Rint_FV2 <- ggplot(data = TableRint, aes(x = Temperature, y = ID_var)) +
  geom_point() + 
  geom_errorbar(aes(ymin = ID_lower, ymax = ID_upper), width = 0) +
  scale_y_continuous(limits = c(0,0.4)) + 
  #scale_x_continuous(breaks = x) + 
  labs(y = "Between-individual variance", x = "Temperature") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

#Plotting FV forest plots - Within ID - amoung series variance in intercept
Rint_FV3 <- ggplot(data = TableRint, aes(x = Temperature, y = Series_var)) +
  geom_point() + 
  geom_errorbar(aes(ymin = Series_lower, ymax = Series_upper), width = 0) +
  scale_y_continuous(limits = c(0,0.4)) + 
  #scale_x_continuous(breaks = x) + 
  labs(y = "Within-individual-between-series\n variance", x = "Temperature") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

#pdf("output/fig/FVvariancecomponentts.pdf", 12, 7)
multiplot(Rint_FV1, Rint_FV2, Rint_FV3, cols = 3)
#dev.off()

#Plotting FV forest plots - Repeatability in Slope
Rslope_FV1 <- ggplot(data = TableRslope, aes(x = Temperature, y = R_slope)) +
  geom_point() + 
  geom_errorbar(aes(ymin = R_lower, ymax = R_upper), width = 0) +
  scale_y_continuous(limits = c(0,1)) + 
  #scale_x_continuous(breaks = x) + 
  labs(y = "Repeatability", x = "Temperature") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

#Plotting FV forest plots - Between ID variance in intercept
Rslope_FV2 <- ggplot(data = TableRslope, aes(x = Temperature, y = ID_var)) +
  geom_point() + 
  geom_errorbar(aes(ymin = ID_lower, ymax = ID_upper), width = 0) +
  scale_y_continuous(limits = c(0,0.2)) + 
  #scale_x_continuous(breaks = x) + 
  labs(y = "Between-individual variance", x = "Temperature") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

#Plotting FV forest plots - Within ID - amoung series variance in intercept
Rslope_FV3 <- ggplot(data = TableRslope, aes(x = Temperature, y = Series_var)) +
  geom_point() + 
  geom_errorbar(aes(ymin = Series_lower, ymax = Series_upper), width = 0) +
  scale_y_continuous(limits = c(0,0.2)) + 
  #scale_x_continuous(breaks = x) + 
  labs(y = "Within-individual-between-series\n variance", x = "Temperature") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

#pdf("output/fig/FV_Rslope_variancecomponentts.pdf", 12, 7)
multiplot(Rslope_FV1, Rslope_FV2, Rslope_FV3, cols = 3)
#dev.off()
