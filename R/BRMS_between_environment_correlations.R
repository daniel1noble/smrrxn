#Using package brms as a sensitivity analysis

setwd("~/Dropbox/smrrxn/")

#clear envir
rm(list = ls())

library(brms)

data <- read.csv("data/data_final/multiresp_log.csv")
data$id <- as.factor(data$id)
data$series <- as.factor(data$series)

complete.data <- data[complete.cases(data),]

fit.t22 <- bf(t_22 ~ z.log.mass + samp_period + (1|p|id))
fit.t24 <- bf(t_24 ~ z.log.mass + samp_period + (1|p|id))
fit.t26 <- bf(t_26 ~ z.log.mass + samp_period + (1|p|id))
fit.t28 <- bf(t_28 ~ z.log.mass + samp_period + (1|p|id))
fit.t30 <- bf(t_30 ~ z.log.mass + samp_period + (1|p|id))
fit.t32 <- bf(t_32 ~ z.log.mass + samp_period + (1|p|id))

brms.m4.usall <- brm(formula = fit.t22 + fit.t24 + fit.t26 + fit.t28 + fit.t30 + fit.t32,
           data = complete.data, family = gaussian(), chains = 2, cores = 2, rescor = T)

brm(fit.t22 +  fit.t24, data = complete.data, family = gaussian())

