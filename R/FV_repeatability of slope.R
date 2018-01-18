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

#m1.log.no.prior
m1.log.no.prior.VCV <- VCV.unpack("output/rds/m1.log.no.prior")
posterior.mode(m1.log.no.prior.VCV)
HPDinterval(as.mcmc(m1.log.no.prior.VCV))

m1.log.no.prior.Sol <- Sol.unpack("output/rds/m1.log.no.prior")
posterior.mode(m1.log.no.prior.Sol)
HPDinterval(as.mcmc(m1.log.no.prior.Sol))

#Repeatability of slope

my.id.slope <- c("inverseK_incb_temp:inverseK_incb_temp.id")
my.series.slope <- c("inverseK_incb_temp:inverseK_incb_temp.series")

rpt.Slope(m1.log.no.prior.VCV, "inverseK_incb_temp:inverseK_incb_temp.id", "inverseK_incb_temp:inverseK_incb_temp.series")





