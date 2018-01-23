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
           data = complete.data, family = gaussian(), chains = 3, cores = 2, iter = 4000, save_model = "stan")
 add_ic(brms.m4.usall) <- c("loo", "waic", "kfold", "R2")
 bayes_R2(brms.m4.usall)
summary(brms.m4.usall)

saveRDS(brms.m4.usall, file = "brms.m4.usall")

## Below some testing to diagnose the problem. Whether the model ran or not is specified above. Problem seems to be specification of rescor argument.
# Model runs
	fit2 <- brm(fit.t22 +  fit.t24, data = complete.data, family = gaussian())

# Model runs
	fit3 <- brm(fit.t22 +  fit.t24, data = complete.data, family = gaussian(), chains = 2, cores = 2)
# Model runs
	fit4 <- brm(fit.t22 +  fit.t24 + fit.t26, data = complete.data, family = gaussian(), chains = 2, cores = 2)
	summary(fit4)

# Model runs
	fit5 <- brm(fit.t22 +  fit.t24 + fit.t26 + fit.t28, data = complete.data, family = gaussian(), chains = 2, cores = 2)
	summary(fit5)

# Model runs
	fit6 <- brm(fit.t22 +  fit.t24 + fit.t26 + fit.t28 + fit.t30, data = complete.data, family = gaussian(), chains = 2, cores = 2)
	summary(fit6)

# Model runs
	fit7 <- brm(fit.t22 + fit.t24 + fit.t26 + fit.t28 + fit.t30 + fit.t32, data = complete.data, family = gaussian(), chains = 2, cores = 2)
	summary(fit7)

# Model doesnt run
	brms.m4.usall <- brm(fit.t22 + fit.t24 + fit.t26 + fit.t28 + fit.t30 + fit.t32,
           data = complete.data, family = gaussian(), chains = 2, cores = 2, rescor = TRUE, save_model = "stan")
# Model runs!
	brms.m4.usall <- brm(formula = fit.t22 + fit.t24 + fit.t26 + fit.t28 + fit.t30 + fit.t32,
           data = complete.data, family = gaussian(), chains = 2, cores = 2, save_model = "stan")