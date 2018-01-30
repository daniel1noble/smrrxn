#Script for create 6 different datasets for releveling the intercept at 6 different temps

setwd("~/Dropbox/smrrxn/")

#clear envir
rm(list = ls())

#libraries
library(dplyr)
library(lme4)

#read in data
data <- read.csv("data/data_final/mrrxn_logT.csv")
data$id <- as.factor(data$id)
data$series <- as.factor(data$series)

varibs.need <- c("samp_period", "id" , "series", "incb_temp","log.temp", "log.prior_temp1",  "log.prior_temp2", "z.log.mass","log.co2pmin", "z.log.co2pmin")

incl.vars <- names(data) %in% varibs.need
data <- data[incl.vars]

str(data)

#The temperatures
sort(unique(data$incb_temp))
sort(unique(data$log.temp))

#Need to center on log level
#current log models will be centered on:
exp(0) #i.e. 1 degrees C

#We want to center them for each model by these values:
sort(unique(data$log.temp))

#Checking
exp(sort(unique(data$log.temp)))


#Centering on log scale
#can't just substract value from 22 on log scale because temps in inverseK scale is not evenly spaced apart. i.e units between 22 - 24 does not equal to 24 -26
tempdat <- data.frame(temp_C = sort(unique(data$incb_temp)),
                      temp_logC = sort(unique(data$log.temp)))

tempdat$dist_log22_to_0 <- c(3.091042 - 3.091042,
                             3.178054 -3.091042,
                             3.258097 - 3.091042,
                             3.332205 - 3.091042,
                             3.401197 -3.091042,
                             3.465736 -3.091042)

tempdat$dist_log24_to_0 <- c(3.091042 - 3.178054,
                             3.178054 -3.178054,
                             3.258097 - 3.178054,
                             3.332205 - 3.178054,
                             3.401197 -3.178054,
                             3.465736 -3.178054)

tempdat$dist_log26_to_0 <- c(3.091042 - 3.258097,
                             3.178054 -3.258097,
                             3.258097 - 3.258097,
                             3.332205 - 3.258097,
                             3.401197 -3.258097,
                             3.465736 -3.258097)

tempdat$dist_log28_to_0 <- c(3.091042 - 3.332205,
                             3.178054 -3.332205,
                             3.258097 - 3.332205,
                             3.332205 - 3.332205,
                             3.401197 -3.332205,
                             3.465736 -3.332205)

tempdat$dist_log30_to_0 <- c(3.091042 - 3.401197,
                             3.178054 -3.401197,
                             3.258097 - 3.401197,
                             3.332205 - 3.401197,
                             3.401197 -3.401197,
                             3.465736 -3.401197)

tempdat$dist_log32_to_0 <- c(3.091042 - 3.465736,
                             3.178054 -3.465736,
                             3.258097 - 3.465736,
                             3.332205 - 3.465736,
                             3.401197 -3.465736,
                             3.465736 -3.465736)

#trying to create variable which specifies 22 is 0
select(tempdat, temp_C, dist_log22_to_0)

vec <- NULL
for(i in 1:length(data$incb_temp)){
  if(data$incb_temp[i] == 22){
    vec[i] <- 0
  } else if(data$incb_temp[i] == 24){
    vec[i] <- 0.087012
  } else if(data$incb_temp[i] == 26){
    vec[i] <- 0.167055
  } else if(data$incb_temp[i] == 28){
    vec[i] <- -0.241163
  } else if(data$incb_temp[i] == 30){
    vec[i] <- 0.310155
  } else if(data$incb_temp[i] == 32){
    vec[i] <- 0.374694
  }
}

str(data) #2520 obs
length(vec) #2520 obs
cbind(vec, data$incb_temp)
data$log.temp_22cen <- vec

#trying to create variable which specifies 24 is 0
select(tempdat, temp_C, dist_log24_to_0)

vec <- NULL
for(i in 1:length(data$incb_temp)){
  if(data$incb_temp[i] == 24){
    vec[i] <- 0
  } else if(data$incb_temp[i] == 22){
    vec[i] <- -0.087012
  } else if(data$incb_temp[i] == 26){
    vec[i] <- 0.080043
  } else if(data$incb_temp[i] == 28){
    vec[i] <- 0.154151
  } else if(data$incb_temp[i] == 30){
    vec[i] <- 0.223143
  } else if(data$incb_temp[i] == 32){
    vec[i] <- 0.287682
  }
}

length(vec) #2520 obs
cbind(vec, data$incb_temp)
data$log.temp_24cen <- vec

#trying to create variable which specifies 26 is 0
select(tempdat, temp_C, dist_log26_to_0)

vec <- NULL
for(i in 1:length(data$incb_temp)){
  if(data$incb_temp[i] == 26){
    vec[i] <- 0
  } else if(data$incb_temp[i] == 22){
    vec[i] <-  -0.167055
  } else if(data$incb_temp[i] == 24){
    vec[i] <- -0.080043
  } else if(data$incb_temp[i] == 28){
    vec[i] <-  0.074108
  } else if(data$incb_temp[i] == 30){
    vec[i] <-  0.143100
  } else if(data$incb_temp[i] == 32){
    vec[i] <- 0.207639
  }
}

length(vec)
cbind(vec, data$incb_temp)
data$log.temp_26cen <- vec

#trying to create variable which specifies 28 is 0
select(tempdat, temp_C, dist_log28_to_0)

vec <- NULL
for(i in 1:length(data$incb_temp)){
  if(data$incb_temp[i] == 28){
    vec[i] <- 0
  } else if(data$incb_temp[i] == 22){
    vec[i] <-  -0.241163
  } else if(data$incb_temp[i] == 24){
    vec[i] <- -0.154151
  } else if(data$incb_temp[i] == 26){
    vec[i] <-  -0.074108
  } else if(data$incb_temp[i] == 30){
    vec[i] <-  0.068992
  } else if(data$incb_temp[i] == 32){
    vec[i] <- 0.133531
  }
}

length(vec)
cbind(vec, data$incb_temp)
data$log.temp_28cen <- vec

#trying to create variable which specifies 30 is 0
select(tempdat, temp_C, dist_log30_to_0)

vec <- NULL
for(i in 1:length(data$incb_temp)){
  if(data$incb_temp[i] == 30){
    vec[i] <- 0
  } else if(data$incb_temp[i] == 22){
    vec[i] <-  -0.310155
  } else if(data$incb_temp[i] == 24){
    vec[i] <- -0.223143
  } else if(data$incb_temp[i] == 26){
    vec[i] <-  -0.143100
  } else if(data$incb_temp[i] == 28){
    vec[i] <- -0.068992
  } else if(data$incb_temp[i] == 32){
    vec[i] <- 0.064539
  }
}

length(vec)
cbind(vec, data$incb_temp)
data$log.temp_30cen <- vec

#trying to create variable which specifies 32 is 0
select(tempdat, temp_C, dist_log32_to_0)

vec <- NULL
for(i in 1:length(data$incb_temp)){
  if(data$incb_temp[i] == 32){
    vec[i] <- 0
  } else if(data$incb_temp[i] == 22){
    vec[i] <- -0.374694
  } else if(data$incb_temp[i] == 24){
    vec[i] <- -0.287682
  } else if(data$incb_temp[i] == 26){
    vec[i] <-  -0.207639
  } else if(data$incb_temp[i] == 28){
    vec[i] <- -0.133531
  } else if(data$incb_temp[i] == 30){
    vec[i] <- -0.064539
  }
}

length(vec)
cbind(vec, data$incb_temp)
data$log.temp_32cen <- vec

write.csv(data, row.names = F, "data/data_final/mr_final_log.T_recentered.csv")

#lme4 model to test recentering of data
str(data)
select(data, incb_temp, log.temp,
       log.temp_22cen, log.temp_24cen, log.temp_26cen, 
       log.temp_28cen, log.temp_30cen, log.temp_32cen) %>% arrange(incb_temp) %>% filter(incb_temp == 32)

#original
model.1 <- lmer(z.log.co2pmin ~ log.temp + z.log.mass + (1+log.temp|id) + (1+log.temp|series), data = data)
summary(model.1)

#22
model.1.22 <- lmer(z.log.co2pmin ~ log.temp_22cen + z.log.mass + (1+log.temp_22cen|id) + (1+log.temp_22cen|series), data = data)
summary(model.1.22)

#24
model.1.24 <- lmer(z.log.co2pmin ~ log.temp_24cen + z.log.mass + (1+log.temp_24cen|id) + (1+log.temp_24cen|series), data = data)
summary(model.1.24)

#26
model.1.26 <- lmer(z.log.co2pmin ~ log.temp_26cen + z.log.mass + (1+log.temp_26cen|id) + (1+log.temp_26cen|series), data = data)
summary(model.1.26)

#28
model.1.28 <- lmer(z.log.co2pmin ~ log.temp_28cen + z.log.mass + (1+log.temp_28cen|id) + (1+log.temp_28cen|series), data = data)
summary(model.1.28)

#30
model.1.30 <- lmer(z.log.co2pmin ~ log.temp_30cen + z.log.mass + (1+log.temp_30cen|id) + (1+log.temp_30cen|series), data = data)
summary(model.1.30)

#32
model.1.32 <- lmer(z.log.co2pmin ~ log.temp_32cen + z.log.mass + (1+log.temp_32cen|id) + (1+log.temp_32cen|series), data = data)
summary(model.1.32)
