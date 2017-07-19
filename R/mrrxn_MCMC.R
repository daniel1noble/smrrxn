setwd("~/Dropbox/smrrxn/")

#clear envir
rm(list = ls())

#load library
library(MCMCglmm)

#read in data

data <- read.csv("data/data_final/mrrxn_final_v2.csv")
data$id <- as.factor(data$id)
data$series <- as.factor(data$series)

varibs.need <- c("obs", "samp_period", "id" ,"batch", "series", "incb_num", "incb_temp_id", "defecate", "incb_temp", "z.incb_temp", "z.log.temp", "incb_temp_K", "inverseK_incb_temp", "body_temp", "z.body_temp", "z.log.body_temp", "body_temp_K", "inverseK_body_temp" , "z.prior_temp1", "z.log.prior_temp1", "prior_temp1_K", "inverseK_prior_temp1", "prior_temp2_K", "inverseK_prior_temp2", "z.prior_temp2", "z.log.prior_temp2", "orig_lizmass", "lizmass_nocombout", "log.mass", "z.log.mass","orig_co2_pmin", "co2pm_nocombout", "log.co2pmin", "z.log.co2pmin")

incl.vars <- names(data) %in% varibs.need
data <- data[incl.vars]

predictors <- c("z.log.co2pmin", "inverseK_incb_temp", "z.log.mass", "inverseK_prior_temp2", "id", "series")

dat <- data[complete.cases(data[,predictors]),]
str(dat)

#priors
expanded.prior <- list(R = list(V = 1, nu = 0.002),
                       G = list(G1 = list(V = diag(2), nu = 0.002, alpha.V = diag(1000,2,2), alpha.mu = rep(0,2)),
                                G2 = list(V = diag(2), nu = 0.002, alpha.V = diag(1000,2,2), alpha.mu = rep(0,2))))

#IW.prior <- list(R = list(V = 1, nu = 0.002),
#                    G = list(G1 = list(V = diag(2), nu = 0.002),
#                             G2 = list(V = diag(2), nu = 0.002)))

#final model - 1 chain for now
#model.1 <- MCMCglmm(z.log.co2pmin ~ inverseK_incb_temp + z.log.mass + inverseK_prior_temp2,
#                    random = ~us(1+inverseK_incb_temp):id + us(1+inverseK_incb_temp):series,
#                    family = "gaussian",
#                    prior = expanded.prior,
#                    nitt = 5010000,
#                    burnin = 10000,
#                    thin = 5000,
#                    data = dat, 
#                    verbose = T)


#saveRDS(model.1, "output/rds/model.1")
#model.1 <- readRDS("output/rds/model.1")

#plot(model.1$Sol)
#summary(model.1$Sol)
#HPDinterval(model.1$Sol)

#heidel.diag(model.1$Sol)  
#geweke.diag(model.1$Sol)

#autocorr.diag(model.1$Sol)
#autocorr.plot(model.1$Sol)

#plot(model.1$VCV)
#summary(model.1$VCV)
#HPDinterval(model.1$VCV)

#heidel.diag(model.1$VCV)  
#geweke.diag(model.1$VCV)

#autocorr.diag(model.1$VCV)
#autocorr.plot(model.1$VCV)


#final model - 1 chain for now
model.2 <- MCMCglmm(z.log.co2pmin ~ inverseK_incb_temp + z.log.mass + inverseK_prior_temp2,
                    random = ~us(1+inverseK_incb_temp):id + us(1+inverseK_incb_temp):series,
                    family = "gaussian",
                    prior = expanded.prior,
                    nitt = 7510000,
                    burnin = 10000,
                    thin = 5000,
                    data = dat, 
                    verbose = T)


#saveRDS(model.2, "output/rds/model.2")
model.2 <- readRDS("output/rds/model.2")

summary(model.2)

#model diagnostics for fixed efs
plot(model.2$Sol)
heidel.diag(model.2$Sol)  
geweke.diag(model.2$Sol)
autocorr.diag(model.2$Sol)
autocorr.plot(model.2$Sol)

#model diagnostics for rand efs
plot(model.2$VCV)
heidel.diag(model.2$VCV)  
geweke.diag(model.2$VCV)
autocorr.diag(model.2$VCV)
autocorr.plot(model.2$VCV)

#output for fix efs
summary(model.2$Sol)
posterior.mode(model.2$Sol)
HPDinterval(model.2$Sol)

#output for ran efs
summary(model.2$VCV)
posterior.mode(model.2$VCV)
HPDinterval(model.2$VCV)

#Tabulating the model output

Table1 <- data.frame(matrix(nrow = 15 , ncol = 3))
rownames(Table1) <- c("Interpcet", "inverseK_incb_temp", "z.log.mass", "inverseK_prior_temp2",
                      "(Intercept).id", "(Intercept):inverseK_incb_temp.id", "inverseK_incb_temp:inverseK_incb_temp.id",
                      "(Intercept).series", "(Intercept):inverseK_incb_temp.series", "inverseK_incb_temp:inverseK_incb_temp.series",
                      "e", "Rint", "Rslope", "Rshort", "Rlong")
colnames(Table1) <- c("estimate", "lower", "upper")

#Tabulating fixed efs and CIs
Table1[1:4, 1] <- posterior.mode(model.2$Sol)
Table1[1:4, 2:3] <- HPDinterval(model.2$Sol)

#Tabulating ran efs and CIs
Table1[5:11, 1] <- posterior.mode(model.2$VCV)[c(1,2,4,5,6,8,9)]
Table1[5:11, 2:3] <- HPDinterval(model.2$VCV)[c(1,2,4,5,6,8,9),]

#Tabulating repeatabilities and CIs
#Rint
Table1[12,1] <- Table1[5,1] / (Table1[5,1] + Table1[8,1]) #Rint
Table1[12,2:3] <- Table1[5,2:3] / (Table1[5,2:3] + Table1[8,2:3]) #Rin CIs

#Rslope
Table1[13,1] <- Table1[6,1] / (Table1[6,1] + Table1[9,1]) #Rslope
Table1[13,2:3] <- Table1[6,2:3] / (Table1[6,2:3] + Table1[9,2:3]) #Rslope CIs

#Rshort
Table1[14,1] <- (Table1[5,1] + Table1[8,1]) / (Table1[5,1] + Table1[8,1] + Table1[11,1]) #Rshort
Table1[14,2:3] <- (Table1[5,2:3] + Table1[8,2:3]) / (Table1[5,2:3] + Table1[8,2:3] + Table1[11,2:3]) #Rshort CIs

#Rlong
Table1[15,1] <- Table1[5,1] / (Table1[5,1] + Table1[8,1] + Table1[11,1]) #Rlong
Table1[15,2:3] <- Table1[5,2:3] / (Table1[5,2:3] + Table1[8,2:3] + Table1[11,2:3]) #Rlong CIs

write.csv(Table1, "output/table/Table1.csv")

#### Pooling three chains

m1 <- mclapply(1:3, function(i) {
  MCMCglmm(z.log.co2pmin ~ inverseK_incb_temp + z.log.mass + inverseK_prior_temp2,
           random = ~us(1+inverseK_incb_temp):id + us(1+inverseK_incb_temp):series,
           family = "gaussian",
           prior = expanded.prior,
           nitt = 75100,
           burnin = 100,
           thin = 50,
           data = dat, 
           verbose = F)
}, mc.cores = 3)


m1.1 <- lapply(m1, function(m) m$Sol)
m1.2 <- do.call(mcmc.list, m1.1)

gelman.plot(m1.2)
gelman.diag(m1.2)
summary(m1.2)
posterior.mode(m1.2)
HPDinterval(as.mcmc(rbind(m1.1[[1]], m1.1[[2]], m1.1[[3]])))

m1.Sol <- lapply(m1, function(m) m$Sol)
m1.Sol <- do.call(mcmc.list, m1.Sol)

m1.VCV <- lapply(m1, function(m) m1.$VCV)
m1.VCV <- do.call(mcmc.list, m1.VCV)

