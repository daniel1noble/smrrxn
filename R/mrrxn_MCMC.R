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
HPDinterval(as.mcmc(rbind(m1.1[[1]], m1.1[[2]], m1.1[[3]])))

#Tabulating the model output

Table1 <- data.frame(matrix(nrow = 15 , ncol = 3))
rownames(Table1) <- c("Interpcet", "inverseK_incb_temp", "z.log.mass", "inverseK_prior_temp2",
                      "(Intercept).id", "(Intercept):inverseK_incb_temp.id", "inverseK_incb_temp:inverseK_incb_temp.id",
                      "(Intercept).series", "(Intercept):inverseK_incb_temp.series", "inverseK_incb_temp:inverseK_incb_temp.series",
                      "e", "Rint", "Rslope", "Rshort", "Rlong")
colnames(Table1) <- c("estimate", "lower", "upper")

#Tabulating fixed efs and CIs
Table1[1:4, 1] <- posterior.mode(m1.Sol)
Table1[1:4, 2:3] <- HPDinterval(as.mcmc(rbind(m1.Sol[[1]], m1.Sol[[2]], m1.Sol[[3]])))

#Tabulating ran efs and CIs
Table1[5:11, 1] <- posterior.mode(m1.VCV)[c(1,2,4,5,6,8,9)]
Table1[5:11, 2:3] <- HPDinterval(as.mcmc(rbind(m1.VCV[[1]], m1.VCV[[2]], m1.VCV[[3]])))[c(1,2,4,5,6,8,9),]

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

source("R/MCMCchains_function.R")

m1 <- MCMC.chains("output/rds/m1")
str(m1)
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
