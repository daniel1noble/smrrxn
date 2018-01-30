#Script for create 6 different datasets for releveling the intercept at 6 different temps

setwd("~/Dropbox/smrrxn/")

#clear envir
rm(list = ls())

#libraries
library(dplyr)

#read in data
data <- read.csv("data/data_final/mrrxn_final_v2.csv")
data$id <- as.factor(data$id)
data$series <- as.factor(data$series)

varibs.need <- c("obs", "samp_period", "id" ,"batch", "series", "incb_num", "incb_temp_id", "defecate", "incb_temp", "z.incb_temp", "z.log.temp", "incb_temp_K", "inverseK_incb_temp", "body_temp", "z.body_temp", "z.log.body_temp", "body_temp_K", "inverseK_body_temp" , "z.prior_temp1", "z.log.prior_temp1", "prior_temp1_K", "inverseK_prior_temp1", "prior_temp2_K", "inverseK_prior_temp2", "z.prior_temp2", "z.log.prior_temp2", "orig_lizmass", "lizmass_nocombout", "log.mass", "z.log.mass","orig_co2_pmin", "co2pm_nocombout", "log.co2pmin", "z.log.co2pmin")

incl.vars <- names(data) %in% varibs.need
data <- data[incl.vars]

str(data)

#The temperatures
unique(data$incb_temp)
unique(data$inverseK_incb_temp)

#Back calculation from inverse K to C

#incb_temp
#data$incb_temp_K <- data$incb_temp + 273.15
#data$inverseK_incb_temp <- 1 / 8.62e-5 * ((1 / mean(data$incb_temp_K)) - (1 /data$incb_temp_K))
#x = -0.6547592
#y = 22
#x = (1 / 8.62e-5) * ((1 / 300.15) - (1 / (y + 273.15)))
#(x / (1 / 8.62e-5)) = ((1 / 300.15) - (1 / (y + 273.15)))
#(x / (1 / 8.62e-5)) - (1 / 300.15) = -(1 / (y + 273.15))
#((x / (1 / 8.62e-5)) - (1 / 300.15))*-1 = (-1 / (y + 273.15))*-1
#((((x / (1 / 8.62e-5)) - (1 / 300.15))*-1) * 273.15) + (y*0.003388108) = 1
#y = (1 - ((((x / (1 / 8.62e-5)) - (1 / 300.15))*-1) * 273.15)) / (((x / (1 / 8.62e-5)) - (1 / 300.15))*-1)

C_to_inverseK <- function(temp_C){
  inverseK <- (1 / 8.62e-5) * ((1 / 300.15) - (1 / (temp_C + 273.15)))
  return(inverseK)
}

inverseK_to_C <- function(inverse_K){
  temp_C <- (1 - ((((inverse_K / (1 / 8.62e-5)) - (1 / 300.15))*-1) * 273.15)) / (((inverse_K / (1 / 8.62e-5)) - (1 / 300.15))*-1)
  return(temp_C)
}


#So at the moment the intercept of the model is inverse_K = 0 
inverseK_to_C(0) 
#which is 27 degrees C using my function
#Example: i need a dataset for 22 degrees 
C_to_inverseK(22) #inverseK for 22degrees

#center on degrees scale first THIS WON'T WORK BECAUSE MODEL WILL STILL SET INTERCEPT AT 0
data$incb_temp22 <- data$incb_temp - 22 #centering on 22 degrees on degrees scale
unique(C_to_inverseK(data$incb_temp22)) #converting deviations from 22 degrees into inverseK
inverseK_to_C(unique(C_to_inverseK(data$incb_temp22))) #back calculating inverseK to degrees

#centering in inverseK scale - NEED TO DO THIS ON THIS LEVEL
#can't just substract value from 22 on inverse scale because temps in inverseK scale is not evenly spaced apart. i.e units between 22 - 24 does not equal to 24 -26
tempdat <- data.frame(temp_C = inverseK_to_C(sort(unique(data$inverseK_incb_temp))),
                      temp_inK = sort(unique(data$inverseK_incb_temp)))

tempdat$dist_frm22_inK_0 <- c(-0.6547592 - -0.6547592,
                              -0.6547592 - -0.3902114,
                              -0.6547592 - -0.1292009,
                              -0.6547592 - 0.1283428,
                              -0.6547592 - 0.3824882,
                              -0.6547592 - 0.6333022)

tempdat$dist_frm24_inK_0 <- c(-0.3902114 - -0.6547592,
                              -0.3902114 - -0.3902114,
                              -0.3902114 - -0.1292009,
                              -0.3902114 - 0.1283428,
                              -0.3902114 - 0.3824882,
                              -0.3902114 - 0.6333022)

tempdat$dist_frm26_inK_0 <- c(-0.1292009 - -0.6547592,
                              -0.1292009 - -0.3902114,
                              -0.1292009 - -0.1292009,
                              -0.1292009 - 0.1283428,
                              -0.1292009 - 0.3824882,
                              -0.1292009 - 0.6333022)

tempdat$dist_frm28_inK_0 <- c(0.1283428 - -0.6547592,
                              0.1283428 - -0.3902114,
                              0.1283428 - -0.1292009,
                              0.1283428 - 0.1283428,
                              0.1283428 - 0.3824882,
                              0.1283428 - 0.6333022)

tempdat$dist_frm30_inK_0 <- c(0.3824882 - -0.6547592,
                              0.3824882 - -0.3902114,
                              0.3824882 - -0.1292009,
                              0.3824882 - 0.1283428,
                              0.3824882 - 0.3824882,
                              0.3824882 - 0.6333022)

tempdat$dist_frm32_inK_0 <- c(0.6333022 - -0.6547592,
                              0.6333022 - -0.3902114,
                              0.6333022 - -0.1292009,
                              0.6333022 - 0.1283428,
                              0.6333022 - 0.3824882,
                              0.6333022 - 0.6333022)

#trying to create variable which specifies 22 is 0
select(tempdat, temp_C, dist_frm22_inK_0)

vec <- NULL
for(i in 1:length(dat$incb_temp)){
  if(dat$incb_temp[i] == 22){
    vec[i] <- 0
  } else if(dat$incb_temp[i] == 24){
    vec[i] <- -0.2645478
  } else if(dat$incb_temp[i] == 26){
    vec[i] <- -0.5255583
  } else if(dat$incb_temp[i] == 28){
    vec[i] <- -0.7831020
  } else if(dat$incb_temp[i] == 30){
    vec[i] <- -1.0372474
  } else if(dat$incb_temp[i] == 32){
    vec[i] <- -1.2880614
  }
}

str(dat) #2401 obs
length(vec) #2401
cbind(vec, dat$incb_temp)
dat$inverseK_incb_temp_22cen <- vec

#trying to create variable which specifies 24 is 0
select(tempdat, temp_C, dist_frm24_inK_0)

vec <- NULL
for(i in 1:length(dat$incb_temp)){
  if(dat$incb_temp[i] == 24){
    vec[i] <- 0
  } else if(dat$incb_temp[i] == 22){
    vec[i] <- 0.2645478
  } else if(dat$incb_temp[i] == 26){
    vec[i] <- -0.2610105
  } else if(dat$incb_temp[i] == 28){
    vec[i] <- -0.5185542
  } else if(dat$incb_temp[i] == 30){
    vec[i] <- -0.7726996
  } else if(dat$incb_temp[i] == 32){
    vec[i] <- -1.0235136
  }
}

length(vec) #2401
cbind(vec, dat$incb_temp)
dat$inverseK_incb_temp_24cen <- vec

#trying to create variable which specifies 26 is 0
select(tempdat, temp_C, dist_frm26_inK_0)

vec <- NULL
for(i in 1:length(dat$incb_temp)){
  if(dat$incb_temp[i] == 26){
    vec[i] <- 0
  } else if(dat$incb_temp[i] == 22){
    vec[i] <-  0.5255583
  } else if(dat$incb_temp[i] == 24){
    vec[i] <- 0.2610105
  } else if(dat$incb_temp[i] == 28){
    vec[i] <- -0.2575437
  } else if(dat$incb_temp[i] == 30){
    vec[i] <- -0.5116891
  } else if(dat$incb_temp[i] == 32){
    vec[i] <- -0.7625031
  }
}

length(vec)
cbind(vec, dat$incb_temp)
dat$inverseK_incb_temp_26cen <- vec

#trying to create variable which specifies 28 is 0
select(tempdat, temp_C, dist_frm28_inK_0)

vec <- NULL
for(i in 1:length(dat$incb_temp)){
  if(dat$incb_temp[i] == 28){
    vec[i] <- 0
  } else if(dat$incb_temp[i] == 22){
    vec[i] <-  0.7831020
  } else if(dat$incb_temp[i] == 24){
    vec[i] <- 0.5185542
  } else if(dat$incb_temp[i] == 26){
    vec[i] <-  0.2575437
  } else if(dat$incb_temp[i] == 30){
    vec[i] <- -0.2541454
  } else if(dat$incb_temp[i] == 32){
    vec[i] <- -0.5049594
  }
}

length(vec)
cbind(vec, dat$incb_temp)
dat$inverseK_incb_temp_28cen <- vec

#trying to create variable which specifies 30 is 0
select(tempdat, temp_C, dist_frm30_inK_0)

vec <- NULL
for(i in 1:length(dat$incb_temp)){
  if(dat$incb_temp[i] == 30){
    vec[i] <- 0
  } else if(dat$incb_temp[i] == 22){
    vec[i] <-  1.0372474
  } else if(dat$incb_temp[i] == 24){
    vec[i] <- 0.7726996
  } else if(dat$incb_temp[i] == 26){
    vec[i] <-  0.5116891
  } else if(dat$incb_temp[i] == 28){
    vec[i] <- 0.2541454
  } else if(dat$incb_temp[i] == 32){
    vec[i] <- -0.2508140
  }
}

length(vec)
cbind(vec, dat$incb_temp)
dat$inverseK_incb_temp_30cen <- vec

#trying to create variable which specifies 32 is 0
select(tempdat, temp_C, dist_frm32_inK_0)

vec <- NULL
for(i in 1:length(dat$incb_temp)){
  if(dat$incb_temp[i] == 32){
    vec[i] <- 0
  } else if(dat$incb_temp[i] == 22){
    vec[i] <- 1.2880614
  } else if(dat$incb_temp[i] == 24){
    vec[i] <- 1.0235136
  } else if(dat$incb_temp[i] == 26){
    vec[i] <-  0.7625031
  } else if(dat$incb_temp[i] == 28){
    vec[i] <- 0.5049594
  } else if(dat$incb_temp[i] == 30){
    vec[i] <- 0.2508140
  }
}

length(vec)
cbind(vec, dat$incb_temp)
dat$inverseK_incb_temp_32cen <- vec

write.csv(dat, row.names = F, "data/data_final/mr_final_recentered.csv")

#lme4 model to test recentering of data
str(dat)
select(dat, incb_temp, inverseK_incb_temp,
       inverseK_incb_temp_22cen, inverseK_incb_temp_24cen, inverseK_incb_temp_26cen, 
       inverseK_incb_temp_28cen, inverseK_incb_temp_30cen, inverseK_incb_temp_32cen) %>% arrange(incb_temp) %>% filter(incb_temp == 32)

#original
model.1 <- lmer(z.log.co2pmin ~ inverseK_incb_temp + z.log.mass + (1+inverseK_incb_temp|id) + (1+inverseK_incb_temp|series), data = dat)
summary(model.1)

#22
model.1.22 <- lmer(z.log.co2pmin ~ inverseK_incb_temp_22cen + z.log.mass + (1+inverseK_incb_temp_22cen|id) + (1+inverseK_incb_temp_22cen|series), data = dat)
summary(model.1.22)

#24
model.1.24 <- lmer(z.log.co2pmin ~ inverseK_incb_temp_24cen + z.log.mass + (1+inverseK_incb_temp_24cen|id) + (1+inverseK_incb_temp_24cen|series), data = dat)
summary(model.1.24)

#26
model.1.26 <- lmer(z.log.co2pmin ~ inverseK_incb_temp_26cen + z.log.mass + (1+inverseK_incb_temp_26cen|id) + (1+inverseK_incb_temp_26cen|series), data = dat)
summary(model.1.26)

#28
model.1.28 <- lmer(z.log.co2pmin ~ inverseK_incb_temp_28cen + z.log.mass + (1+inverseK_incb_temp_28cen|id) + (1+inverseK_incb_temp_28cen|series), data = dat)
summary(model.1.28)

#30
model.1.30 <- lmer(z.log.co2pmin ~ inverseK_incb_temp_30cen + z.log.mass + (1+inverseK_incb_temp_30cen|id) + (1+inverseK_incb_temp_30cen|series), data = dat)
summary(model.1.30)

#32
model.1.32 <- lmer(z.log.co2pmin ~ inverseK_incb_temp_32cen + z.log.mass + (1+inverseK_incb_temp_32cen|id) + (1+inverseK_incb_temp_32cen|series), data = dat)
summary(model.1.32)
