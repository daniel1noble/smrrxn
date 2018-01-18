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

##Yimen method
#m1.log.np.t22
VCV.check("output/rds/m1.log.np.t22")
Sol.check("output/rds/m1.log.np.t22")

m1.log.np.t22.VCV <- VCV.unpack("output/rds/m1.log.np.t22")
posterior.mode(m1.log.np.t22.VCV)
HPDinterval(as.mcmc(m1.log.np.t22.VCV))

#Repeatability of intercept
FV.rpt.Int(m1.log.np.t22.VCV)

#Repeatability 'short term'
FV.rpt.Int.short(m1.log.np.t22.VCV)

#Repeatability 'long term' - classic repeatability
FV.rpt.Int.long(m1.log.np.t22.VCV)

#m1.log.np.t24
VCV.check("output/rds/m1.log.np.t24")
Sol.check("output/rds/m1.log.np.t24")

m1.log.np.t24.VCV <- VCV.unpack("output/rds/m1.log.np.t24")
colnames(m1.log.np.t24.VCV)[c(1,5)]
posterior.mode(m1.log.np.t24.VCV)
HPDinterval(as.mcmc(m1.log.np.t24.VCV))

#Repeatability of intercept
FV.rpt.Int(m1.log.np.t24.VCV)

#Repeatability 'short term'
FV.rpt.Int.short(m1.log.np.t24.VCV)

#Repeatability 'long term' - classic repeatability
FV.rpt.Int.long(m1.log.np.t24.VCV)

#m1.log.np.t26
VCV.check("output/rds/m1.log.np.t26")
Sol.check("output/rds/m1.log.np.t26")

m1.log.np.t26.VCV <- VCV.unpack("output/rds/m1.log.np.t26")
colnames(m1.log.np.t26.VCV)[c(1,5)]
posterior.mode(m1.log.np.t26.VCV)
HPDinterval(as.mcmc(m1.log.np.t26.VCV))

#Repeatability of intercept
FV.rpt.Int(m1.log.np.t26.VCV)

#Repeatability 'short term'
FV.rpt.Int.short(m1.log.np.t26.VCV)

#Repeatability 'long term' - classic repeatability
FV.rpt.Int.long(m1.log.np.t26.VCV)

#m1.log.np.t28
VCV.check("output/rds/m1.log.np.t28")
Sol.check("output/rds/m1.log.np.t28")

m1.log.np.t28.VCV <- VCV.unpack("output/rds/m1.log.np.t28")
colnames(m1.log.np.t28.VCV)[c(1,5)]
posterior.mode(m1.log.np.t28.VCV)
HPDinterval(as.mcmc(m1.log.np.t28.VCV))

#Repeatability of intercept
FV.rpt.Int(m1.log.np.t28.VCV)

#Repeatability 'short term'
FV.rpt.Int.short(m1.log.np.t28.VCV)

#Repeatability 'long term' - classic repeatability
FV.rpt.Int.long(m1.log.np.t28.VCV)

#m1.log.np.t30
VCV.check("output/rds/m1.log.np.t30")
Sol.check("output/rds/m1.log.np.t30")

m1.log.np.t30.VCV <- VCV.unpack("output/rds/m1.log.np.t30")
colnames(m1.log.np.t30.VCV)[c(1,5)]
posterior.mode(m1.log.np.t30.VCV)
HPDinterval(as.mcmc(m1.log.np.t30.VCV))

#Repeatability of intercept
FV.rpt.Int(m1.log.np.t30.VCV)

#Repeatability 'short term'
FV.rpt.Int.short(m1.log.np.t30.VCV)

#Repeatability 'long term' - classic repeatability
FV.rpt.Int.long(m1.log.np.t30.VCV)

#m1.log.np.t32
VCV.check("output/rds/m1.log.np.t32")
Sol.check("output/rds/m1.log.np.t32")

m1.log.np.t32.VCV <- VCV.unpack("output/rds/m1.log.np.t32")
colnames(m1.log.np.t32.VCV)[c(1,5)]
posterior.mode(m1.log.np.t32.VCV)
HPDinterval(as.mcmc(m1.log.np.t32.VCV))

#Repeatability of intercept
FV.rpt.Int(m1.log.np.t32.VCV)

#Repeatability 'short term'
FV.rpt.Int.short(m1.log.np.t32.VCV)

#Repeatability 'long term' - classic repeatability
FV.rpt.Int.long(m1.log.np.t32.VCV)


##Singer and Willet method
#m1.log.noprior.noseries
VCV.check("output/rds/m1.log.noprior.noseries")
Sol.check("output/rds/m1.log.noprior.noseries")

m1.log.np.ns.VCV <- VCV.unpack("output/rds/m1.log.noprior.noseries")
colnames(m1.log.np.ns.VCV)

#Repeatabilty at 22
SW.rpt.Temp(m1.log.np.ns.VCV, temp = 22)

#Repeatabilty at 24
SW.rpt.Temp(m1.log.np.ns.VCV, temp = 24)

#Repeatabilty at 26
SW.rpt.Temp(m1.log.np.ns.VCV, temp = 26)
  
#Repeatabilty at 28
SW.rpt.Temp(m1.log.np.ns.VCV, temp = 28)  

#Repeatabilty at 30
SW.rpt.Temp(m1.log.np.ns.VCV, temp = 30)

#Repeatabilty at 32
SW.rpt.Temp(m1.log.np.ns.VCV, temp = 32)

##Character-state way
#m4.log.usall

VCV.check("output/rds/m4.log.usall")
Sol.check("output/rds/m4.log.usall")

m4.log.usall.VCV <- VCV.unpack("output/rds/m4.log.usall")

#Repeatability at 22
CS.rpt.Temp(m4.log.usall.VCV, temp = 22)

CS.rpt.Temp(m4.log.usall.VCV, temp = 24)

CS.rpt.Temp(m4.log.usall.VCV, temp = 26)

CS.rpt.Temp(m4.log.usall.VCV, temp = 28)

CS.rpt.Temp(m4.log.usall.VCV, temp = 30)

CS.rpt.Temp(m4.log.usall.VCV, temp = 32)
















