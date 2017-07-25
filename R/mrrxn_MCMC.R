setwd("~/Dropbox/smrrxn/")

#clear envir
rm(list = ls())

#load library
library(MCMCglmm)

m1 <- readRDS("output/rds/m1")

m1.S <- lapply(m1, function(m) m$Sol)
m1.Sol <- do.call(mcmc.list, m1.S)

m1.V <- lapply(m1, function(m) m$VCV)
m1.VCV <- do.call(mcmc.list, m1.V)

gelman.diag(m1.Sol)
summary(m1.Sol)
posterior.mode(m1.Sol)
HPDinterval(as.mcmc(rbind(m1.Sol[[1]], m1.Sol[[2]], m1.Sol[[3]])))

gelman.diag(m1.VCV, multivariate = F)
summary(m1.VCV)
posterior.mode(m1.VCV)
HPDinterval(as.mcmc(rbind(m1.VCV[[1]], m1.VCV[[2]], m1.VCV[[3]])))

#Tabulating the model output

Table1 <- data.frame(matrix(nrow = 15 , ncol = 3))
#rownames(Table1) <- c("Interpcet", "inverseK_incb_temp", "z.log.mass", "inverseK_prior_temp2",
#                      "IDintercept", "IDslope", "COV IDintercept-IDslope",
#                      "Sintercept", "Sslope", "COV Sintercept-Sslope",
#                      "e", "Rint", "Rslope", "Rshort", "Rlong")

rownames(Table1) <- c(names(c(posterior.mode(m1.Sol), posterior.mode(m1.VCV)[c(1,4,2,5,8,6,9)])), "Rint", "Rslope", "Rshort", "Rlong")
colnames(Table1) <- c("estimate", "lower", "upper")

#Tabulating fixed efs and CIs
Table1[1:4, 1] <- posterior.mode(m1.Sol)
Table1[1:4, 2:3] <- HPDinterval(as.mcmc(rbind(m1.Sol[[1]], m1.Sol[[2]], m1.Sol[[3]])))

#Tabulating ran efs and CIs
Table1[5:11, 1] <- posterior.mode(m1.VCV)[c(1,4,2,5,8,6,9)]
Table1[5:11, 2:3] <- HPDinterval(as.mcmc(rbind(m1.VCV[[1]], m1.VCV[[2]], m1.VCV[[3]])))[c(1,4,2,5,8,6,9),]

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
write.csv(round(Table1,2), "output/table/Table1_rounded.csv")

#trying to make predictions

#source("R/MCMCchains_function.R")

#m2 <- MCMC.chains("output/rds/m2")

m2 <- readRDS("output/rds/m2")

m2.S <- lapply(m2, function(m) m$Sol)
m2.Sol <- do.call(mcmc.list, m2.S)

plot(m2.Sol)
gelman.diag(m2.Sol)
summary(m2.Sol)
posterior.mode(m2.Sol)
HPDinterval(as.mcmc(rbind(m2.Sol[[1]], m2.Sol[[2]], m2.Sol[[3]])))

m2.V <- lapply(m2, function(m) m$VCV)
m2.VCV <- do.call(mcmc.list, m2.V)

plot(m2.VCV)
gelman.diag(m2.VCV, multivariate = F)
summary(m2.VCV)
posterior.mode(m2.VCV)
HPDinterval(as.mcmc(rbind(m2.VCV[[1]], m2.VCV[[2]], m2.VCV[[3]])))


str(m2)
str(m2$solVCVlist)
names(data) 

newdata = data.frame(inverseK_incb_temp = rep(c(-0.6547592,  0.6333022, -0.1292009,  0.1283428, -0.3902114,  0.3824882),42*10),
                     z.log.mass = rep(0, 420*6),
                     inverseK_prior_temp2 = rep(-0.001947, 420*6),
                     id = rep(unique(data$id), 10*6),
                     series = rep(seq(1:10), each = 42*6))

X <-model.matrix(~inverseK_incb_temp + z.log.mass + inverseK_prior_temp2 , data=newdata)

beta <- m1$solVCVchain$Sol[,1]

newdata$pred_zco2pm <- X %*% beta
newdata$pred_logco2pm <-(newdata$pred_zco2pm * sd(newdata$pred_zco2pm)) + mean(newdata$pred_zco2pm)
newdata$pred_co2pm <- exp(newdata$pred_logco2pm)


predict.MCMCglmm(m1, newdata = newdata)
