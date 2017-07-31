setwd("~/Dropbox/smrrxn/")

#clear envir
rm(list = ls())

#load library
library(MCMCglmm)

# Define whether to run code

m1 <- FALSE
m2 <- FALSE
m3 <- FALSE

#priors
expanded.prior <- list(R = list(V = 1, nu = 0.002),
                       G = list(G1 = list(V = diag(2), nu = 0.002, alpha.V = diag(1000,2,2), alpha.mu = rep(0,2)),
                                G2 = list(V = diag(2), nu = 0.002, alpha.V = diag(1000,2,2), alpha.mu = rep(0,2))))

#m1
if(m1){
  m1 <- mclapply(1:3, function(i) {
    MCMCglmm(z.log.co2pmin ~ z.log.co2pmin ~ inverseK_incb_temp + z.log.mass + inverseK_prior_temp2,
             random = ~us(1+inverseK_incb_temp):id + us(1+inverseK_incb_temp):series,
             family = "gaussian",
             prior = expanded.prior,
             nitt = 5010000,
             burnin = 10000,
             thin = 5000,
             data = dat, 
             verbose = T)
  }, mc.cores = 3)
}else{
m1 <- readRDS("output/rds/m1")
}

m1.S <- lapply(m1, function(m) m$Sol)
m1.Sol <- do.call(rbind, m1.S)

m1.V <- lapply(m1, function(m) m$VCV)
m1.VCV <- do.call(rbind, m1.V)

gelman.diag(m1.Sol)
summary(m1.Sol)
posterior.mode(m1.Sol)
HPDinterval(as.mcmc(rbind(m1.Sol[[1]], m1.Sol[[2]], m1.Sol[[3]])))

gelman.diag(m1.VCV, multivariate = F)
summary(m1.VCV)
posterior.mode(m1.VCV)
HPDinterval(as.mcmc(rbind(m1.VCV[[1]], m1.VCV[[2]], m1.VCV[[3]])))

#Calculating repeatabilities using m2
#m2
if(m2){
  m2 <- mclapply(1:3, function(i) {
    MCMCglmm(z.log.co2pmin ~ z.log.co2pmin ~ inverseK_incb_temp + z.log.mass + inverseK_prior_temp2,
             random = ~us(1+inverseK_incb_temp):id + us(1+inverseK_incb_temp):series,
             family = "gaussian",
             prior = expanded.prior,
             nitt = 7510000,
             burnin = 10000,
             thin = 5000,
             data = dat, 
             verbose = T)
  }, mc.cores = 3)
}else{
  m2 <- readRDS("output/rds/m2")
}

m2.S <- lapply(m2, function(m) m$Sol)
m2.S1 <- do.call(rbind, m2.S)

#Tabulating the model output using m2
Table1 <- data.frame(matrix(nrow = 15 , ncol = 3))
#rownames(Table1) <- c("Interpcet", "inverseK_incb_temp", "z.log.mass", "inverseK_prior_temp2",
#                      "IDintercept", "IDslope", "COV IDintercept-IDslope",
#                      "Sintercept", "Sslope", "COV Sintercept-Sslope",
#                      "e", "Rint", "Rslope", "Rshort", "Rlong")

rownames(Table1) <- c(names(c(posterior.mode(m2.Sol), posterior.mode(m2.VCV)[c(1,4,2,5,8,6,9)])), "Rint", "Rslope", "Rshort", "Rlong")
colnames(Table1) <- c("estimate", "lower", "upper")

#Tabulating fixed efs and CIs
Table1[1:4, 1] <- posterior.mode(m2.Sol)
Table1[1:4, 2:3] <- HPDinterval(as.mcmc(rbind(m2.Sol[[1]], m2.Sol[[2]], m2.Sol[[3]])))

#Tabulating ran efs and CIs
Table1[5:11, 1] <- posterior.mode(m2.VCV)[c(1,4,2,5,8,6,9)]
Table1[5:11, 2:3] <- HPDinterval(as.mcmc(rbind(m2.VCV[[1]], m2.VCV[[2]], m2.VCV[[3]])))[c(1,4,2,5,8,6,9),]

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

#Individual predictions from m2 for plots

post <- m2.S1

extractID <- function(post, fixef = colnames(post)[1:4], id = "ld0133"){
    columns <- as.character(c(fixef, colnames(post)[grep(id, colnames(post))]))  
    solID <- post[,columns]
    return(solID)
}
test <- extractID(m2.S1)

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

#Intercept model for calculating PCV 
#m3
if(m3){
  m3 <- mclapply(1:3, function(i) {
    MCMCglmm(z.log.co2pmin ~ 1,
             random = ~us(1+inverseK_incb_temp):id + us(1+inverseK_incb_temp):series,
             family = "gaussian",
             prior = expanded.prior,
             nitt = 7510000,
             burnin = 10000,
             thin = 5000,
             data = dat, 
             verbose = T)
  }, mc.cores = 3)
}else{
  m3 <- readRDS("output/rds/m3")
}

