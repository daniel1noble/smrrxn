setwd("~/Dropbox/smrrxn/")

#clear envir
rm(list = ls())

#load library
library(MCMCglmm)
library(parallel)
library(dplyr)
library(ggplot2)
library(reshape2)

source("R/functions/smr_functions.R")
source("R/functions/unpacking.R")

#m1.log
VCV.check("output/rds/m1.log")
m1.log.VCV <- VCV.unpack("output/rds/m1.log")
posterior.mode(m1.log.VCV)
HPDinterval(as.mcmc(m1.log.VCV))

Sol.check("output/rds/m1.log")
m1.log.Sol <- Sol.unpack("output/rds/m1.log")
posterior.mode(m1.log.Sol)
HPDinterval(as.mcmc(m1.log.Sol))

#DIC - #m1.log 1657.663
mean(unlist(my.DIC("output/rds/m1.log")))

#R^2 - #m1.log
m1.log <- readRDS("output/rds/m1.log")

m1.log.mFixed <- mean(m1.log[[1]]$Sol[,2]) * m1.log[[1]]$X[,2] + mean(m1.log[[1]]$Sol[,3]) * m1.log[[1]]$X[,3] + mean(m1.log[[1]]$Sol[,4]) * m1.log[[1]]$X[,4] 

m1.log.mVarF <- var(m1.log.mFixed)
m1.log.mVarF <- var(as.vector(apply(m1.log[[1]]$Sol,2,mean) %*% t(m1.log[[1]]$X)))

#m1.log MCMCglmm - marginal 0.4304141
m1.log.mVarF/(m1.log.mVarF + sum(apply(m1.log[[1]]$VCV,2,mean)))

#m1.log.no.prior
VCV.check("output/rds/m1.log.no.prior")
m1.log.no.prior.VCV <- VCV.unpack("output/rds/m1.log.no.prior")
posterior.mode(m1.log.no.prior.VCV)
HPDinterval(as.mcmc(m1.log.no.prior.VCV))

Sol.check("output/rds/m1.log.no.prior")
m1.log.no.prior.Sol <- Sol.unpack("output/rds/m1.log.no.prior")
posterior.mode(m1.log.no.prior.Sol)
HPDinterval(as.mcmc(m1.log.no.prior.Sol))

#DIC - m1.log.no.prior 1656.23
mean(unlist(my.DIC("output/rds/m1.log.no.prior")))

#R^2 - m1.log.no.prior
m1.log.no.prior <- readRDS("output/rds/m1.log.no.prior")

m1.log.no.prior.mFixed <- mean(m1.log.no.prior[[1]]$Sol[,2]) * m1.log.no.prior[[1]]$X[,2] + mean(m1.log.no.prior[[1]]$Sol[,3]) * m1.log.no.prior[[1]]$X[,3]

m1.log.no.prior.mVarF <- var(m1.log.no.prior.mFixed)
m1.log.no.prior.mVarF <- var(as.vector(apply(m1.log.no.prior[[1]]$Sol,2,mean) %*% t(m1.log.no.prior[[1]]$X)))

#m1.log.no.prior MCMCglmm - marginal 0.4287207
m1.log.no.prior.mVarF/(m1.log.no.prior.mVarF + sum(apply(m1.log.no.prior[[1]]$VCV,2,mean)))
