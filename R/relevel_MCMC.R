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

#Plotting these temp relationships out
unique(data$incb_temp)
unique(data$inverseK_incb_temp)

par(mfrow = c(1,2))

mod.inK <- lm(z.log.co2pmin ~ inverseK_incb_temp, data = data)
mod.C <- lm(z.log.co2pmin ~ incb_temp, data = data)

plot(data$incb_temp, data$z.log.co2pmin)
abline(mod.C)
plot(data$inverseK_incb_temp, data$z.log.co2pmin)
abline(mod.inK)


#lme4 model to test recentering of data
model.1 <- lmer(z.log.co2pmin ~ inverseK_incb_temp + z.log.mass + (1+inverseK_incb_temp|id) + (1+inverseK_incb_temp|series), data = data)
summary(model.1)

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

filter(dat, incb_temp == 22) %>% nrow +
filter(dat, incb_temp == 24) %>% nrow +
filter(dat, incb_temp == 26) %>% nrow +
filter(dat, incb_temp == 28) %>% nrow +
filter(dat, incb_temp == 30) %>% nrow +
filter(dat, incb_temp == 32) %>% nrow
  
length(vec)

head(vec) 
head(dat$incb_temp)

cbind(vec, dat$incb_temp)
