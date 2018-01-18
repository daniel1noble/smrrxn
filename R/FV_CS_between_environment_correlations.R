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

##Character-state way
#m4.log.usall

VCV.check("output/rds/m4.log.usall")
Sol.check("output/rds/m4.log.usall")

m4.log.usall.VCV <- VCV.unpack("output/rds/m4.log.usall")

#between covariances and correlations
#Cov not positive definitive to convert to correlation
my.cor.cov.matrices(m4.log.usall.VCV, type = "between")

#within covariances and correlations
my.cor.cov.matrices(m4.log.usall.VCV, type = "within")

#Losing power, so will estimate covariation at between level only
#m4.log.idhwith

VCV.check("output/rds/m4.log.idhwith")
Sol.check("output/rds/m4.log.idhwith")

m4.log.idhwith.VCV <- VCV.unpack("output/rds/m4.log.idhwith")

#between covariances and correlations
#Cov not positive definitive, didn't convert to correlation correctly
my.cor.cov.matrices(m4.log.idhwith.VCV, type = "between")

#Lets try this interaction fixed effect model with idhwith
#m4.interaction.idhwith

VCV.check("output/rds/m4.interaction.idhwith")
Sol.check("output/rds/m4.interaction.idhwith")

m4.interaction.idhwith.VCV <- VCV.unpack("output/rds/m4.interaction.idhwith")

#between covariances and correlations
#Oddly, estimating more fixed effects...but the covariance matrix is okay to convert and the direction of covariation has opposite relative to previous two models
my.cor.cov.matrices(m4.interaction.idhwith.VCV, type = "between")

#Lets try estimate everything include us(with)
#m4.interaction.usall

VCV.check("output/rds/m4.interaction.usall")
Sol.check("output/rds/m4.interaction.usall")

m4.interaction.usall.VCV <- VCV.unpack("output/rds/m4.interaction.usall")

#between covariances and correlations - no issues with converting Cov to Cor
my.cor.cov.matrices(m4.interaction.usall.VCV, type = "between")
