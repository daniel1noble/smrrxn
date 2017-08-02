setwd("~/Dropbox/smrrxn/")

#clear envir
rm(list = ls())

#load library
library(MCMCglmm)
library(parallel)

#data

data <- read.csv("data/data_final/mrrxn_final_v2.csv")
data$id <- as.factor(data$id)
data$series <- as.factor(data$series)

varibs.need <- c("obs", "samp_period", "id" ,"batch", "series", "incb_num", "incb_temp_id", "defecate", "incb_temp", "z.incb_temp", "z.log.temp", "incb_temp_K", "inverseK_incb_temp", "body_temp", "z.body_temp", "z.log.body_temp", "body_temp_K", "inverseK_body_temp" , "z.prior_temp1", "z.log.prior_temp1", "prior_temp1_K", "inverseK_prior_temp1", "prior_temp2_K", "inverseK_prior_temp2", "z.prior_temp2", "z.log.prior_temp2", "orig_lizmass", "lizmass_nocombout", "log.mass", "z.log.mass","orig_co2_pmin", "co2pm_nocombout", "log.co2pmin", "z.log.co2pmin")

incl.vars <- names(data) %in% varibs.need
data <- data[incl.vars]

predictors <- c("z.log.co2pmin", "inverseK_incb_temp", "z.log.mass", "inverseK_prior_temp2", "id", "series")

dat <- data[complete.cases(data[,predictors]),]
str(dat)

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

#Calculating repeatabilities using m1
m1.S <- lapply(m1, function(m) m$Sol)
m1.Sol <- do.call(rbind, m1.S)
m1.S1 <- do.call(mcmc.list, m1.S)

m1.V <- lapply(m1, function(m) m$VCV)
m1.VCV <- do.call(rbind, m1.V)
m1.V1 <- do.call(mcmc.list, m1.V)

gelman.diag(m1.S1)
summary(m1.Sol)
posterior.mode(m1.Sol)
HPDinterval(as.mcmc(m1.Sol))

gelman.diag(m1.V1, multivariate = F)
summary(m1.VCV)
posterior.mode(m1.VCV)
HPDinterval(as.mcmc(m1.VCV))

#Tabulating the model output using m1
Table1 <- data.frame(matrix(nrow = 15 , ncol = 3))
#rownames(Table1) <- c("Interpcet", "inverseK_incb_temp", "z.log.mass", "inverseK_prior_temp2",
#                      "IDintercept", "IDslope", "COV IDintercept-IDslope",
#                      "Sintercept", "Sslope", "COV Sintercept-Sslope",
#                      "e", "Rint", "Rslope", "Rshort", "Rlong")

rownames(Table1) <- c(names(c(posterior.mode(m1.Sol), posterior.mode(m1.VCV)[c(1,4,2,5,8,6,9)])), "Rint", "Rslope", "Rshort", "Rlong")
colnames(Table1) <- c("estimate", "lower", "upper")

#Tabulating fixed efs and CIs
Table1[1:4, 1] <- posterior.mode(m1.Sol)
Table1[1:4, 2:3] <- HPDinterval(as.mcmc(m1.Sol))

#Tabulating ran efs and CIs
Table1[5:11, 1] <- posterior.mode(m1.VCV)[c(1,4,2,5,8,6,9)]
Table1[5:11, 2:3] <- HPDinterval(as.mcmc(m1.VCV))[c(1,4,2,5,8,6,9),]

#Tabulating repeatabilities and CIs
head(m1.VCV)
colnames(m1.VCV)

#Rint
R.int <- m1.VCV[,"(Intercept):(Intercept).id"] / ( m1.VCV[,"(Intercept):(Intercept).id"] + m1.VCV[,"(Intercept):(Intercept).series"] ) 

Table1[12,1] <- posterior.mode(R.int) #R.int
Table1[12,2:3] <- HPDinterval(as.mcmc(R.int)) #R.int CIs

#Rslope
R.slope <- m1.VCV[,"inverseK_incb_temp:inverseK_incb_temp.id"] / ( m1.VCV[,"inverseK_incb_temp:inverseK_incb_temp.id"] + m1.VCV[,"inverseK_incb_temp:inverseK_incb_temp.series"] )

Table1[13,1] <-  posterior.mode(R.slope) #Rslope
Table1[13,2:3] <- HPDinterval(as.mcmc(R.slope)) #Rslope CIs

#Rshort
R.short <- ( m1.VCV[,"(Intercept):(Intercept).id"] + m1.VCV[,"(Intercept):(Intercept).series"] ) / ( m1.VCV[,"(Intercept):(Intercept).id"] + m1.VCV[,"(Intercept):(Intercept).series"]  + m1.VCV[,"units"] )

Table1[14,1] <- posterior.mode(R.short) #R.short
Table1[14,2:3] <- HPDinterval(as.mcmc(R.short)) #R.short CIs

#Rlong
R.long <-  m1.VCV[,"(Intercept):(Intercept).id"] / ( m1.VCV[,"(Intercept):(Intercept).id"] + m1.VCV[,"(Intercept):(Intercept).series"]  + m1.VCV[,"units"] )

Table1[15,1] <-  posterior.mode(R.long) #Rlong
Table1[15,2:3] <- HPDinterval(as.mcmc(R.long)) #Rlong CIs

write.csv(Table1, "output/table/Table1.csv")
write.csv(round(Table1,2), "output/table/Table1_rounded.csv")

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
m2.Sol <- do.call(rbind, m2.S)
m2.S1 <- do.call(mcmc.list, m2.S)

gelman.diag(m2.S1, multivariate = F)
summary(m2.Sol)
posterior.mode(m2.Sol)
HPDinterval(as.mcmc(m2.Sol))

m2.V <- lapply(m2, function(m) m$VCV)
m2.V1 <- do.call(mcmc.list, m2.V)
m2.VCV <- do.call(rbind, m2.V)

plot(m2.VCV)
gelman.diag(m2.VCV, multivariate = F)
summary(m2.VCV)
posterior.mode(m2.VCV)
HPDinterval(as.mcmc(rbind(m2.VCV[[1]], m2.VCV[[2]], m2.VCV[[3]])))


#Individual predictions from m2 for plots

post <- m2.S1

extractID <- function(post, fixef = colnames(post)[1:4], id = "ld0133"){
    columns <- as.character(c(fixef, colnames(post)[grep(id, colnames(post))]))  
    solID <- post[,columns]
    return(solID)
}
test <- extractID(m2.S1)

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

