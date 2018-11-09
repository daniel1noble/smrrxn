setwd("~/gitrepo/smrrxn/")

#clear envir
rm(list = ls())

#load library
library(MCMCglmm)
library(parallel)

#data
data <- read.csv("data/data_final/multiresp_log_pT.csv")
data$id <- as.factor(data$id)
data$series <- as.factor(data$series)

data <- read.csv("data/data_final/mrrxn_logT.csv")

data3 <- data %>% 
select(id, samp_period, series, incb_temp, log.co2pmin, orig_lizmass, prior_temp1) %>%
group_by(series) %>%
mutate(avg_mass = mean(orig_lizmass), avg_log_prior_temp1 = log(mean(prior_temp1, na.rm = T))) %>%
select(id, samp_period, series, avg_mass, avg_log_prior_temp1, incb_temp, log.co2pmin) %>%
spread(incb_temp, log.co2pmin) 


colnames(data3)[6:11] <- paste0("t_", colnames(data3)[6:11])

data3$z.log.mass <- scale(log(data3$avg_mass))

data3 <- as.data.frame(data3)
data3 <- data3[!is.na(data3$z.log.mass),]
write.csv(data3, row.names = F, "data/data_final/multiresp_log_pT.csv")

#priors
multi.prior <- list(R = list(V = diag(6), nu = 0.01), 
                    G = list(G1 = list(V = diag(6), nu = 0.01),
                             G2 = list(V = diag(6), nu = 0.01)))

#m4.logMR.idhsesh.usall
m4.logMR.idhsesh.usall <- mclapply(1:3, function(i) {
  MCMCglmm(cbind(t_22, t_24, t_26, t_28, t_30, t_32) ~  z.log.mass + avg_log_prior_temp1,
           random= ~us(trait):id + idh(trait):samp_period,
           rcov = ~us(trait):units,
           family = c(rep("gaussian", 6)),
           prior = multi.prior,
           nitt = 7510,
           burnin = 100,
           thin = 50,
           data = data, 
           verbose = T)
}, mc.cores = 3)

saveRDS(m4.logMR.idhsesh.usall, "R/m4.logMR.idhsesh.usall")






