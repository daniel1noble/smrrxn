setwd("~/Dropbox/smrrxn/")

rm(list = ls())

library(ggplot2)
library(lme4)

rawdata <- read.csv("data/data_final/smr_rawfinal.csv")
str(rawdata)

rawdata <- rawdata[order(rawdata$samp_period),]

data <- rawdata[,-c(12,13,15,16,17,18,19,20)]
names(data)
str(data)

#plotting to see if there is non-linear relationship
boxplot(data$co2_pmin ~ data$incb_temp) #Doesn't look like it

plot(log(data$co2_pmin) ~ log(data$lizmass)) #Doesn't look like it

#Outlier in body temp
ggplot(data, aes(body_temp, obs, label = obs)) +
  geom_point(colour = "white") +
  geom_text(check_overlap = T) + 
  theme_bw()

data[data$obs == 1394,12] <- 24
data[data$obs == 1399,12] <- 24 #No more outliers in body temp

#Collinearity between temp and mass? 

pairs(data[,c(25,10,22)], log = "xy")
cor(data[,c(25,10,22,11)],use = "complete.obs") #thats good no super highly correlated variables

#standardising variables

data$z.incb_temp <- scale(data$incb_temp)
data$z.body_temp <- scale(data$body_temp)

data$log.mass <- log(data$lizmass)
data$z.log.mass <- scale(data$log.mass)

data$log.co2pmin <- log(data$co2_pmin)
data$z.log.co2pmin <- scale(data$log.co2pmin)

write.csv(data, "data/data_final/mr_data_analysis.csv", row.names = F)

#histograms

data <- read.csv("data/data_final/mr_data_analysis_example.csv")

hist(data$z.incb_temp) ; shapiro.test(data$z.incb_temp)
hist(log(data$incb_temp)) ; shapiro.test(log(data$incb_temp))

hist(data$z.body_temp) ; shapiro.test(data$z.body_temp)

hist(data$log.mass) ; shapiro.test(data$log.mass)
hist(data$z.log.mass) ; shapiro.test(data$z.log.mass)

hist(data$log.co2pmin) ; shapiro.test(data$log.co2pmin)
hist(data$z.log.co2pmin) ; shapiro.test(data$z.log.co2pmin)

#batch and incubation number and defecation effects
#on mr
boxplot(data$z.log.co2pmin ~ data$site)
boxplot(data$z.log.co2pmin ~ data$batch)
boxplot(data$z.log.co2pmin ~ data$incb_num)
boxplot(data$z.log.co2pmin ~ data$defecate)

#on mass
boxplot(data$z.log.mass ~ data$site) #Quite diff in mass between pops
boxplot(data$z.log.mass ~ data$batch)
boxplot(data$z.log.mass ~ data$incb_num)
boxplot(data$z.log.mass ~ data$defecate)

#Carry over effects of previous temperature
data$z.prior_temp_1 <- scale(data$prior_temp_1)
data$z.prior_temp_2 <- scale(data$prior_temp_2)
data$z.prior_temp_3 <- scale(data$prior_temp_3)
data$z.prior_temp_4 <- scale(data$prior_temp_4)

data$z.log.prior_temp_1 <- scale(log(data$prior_temp_1))
data$z.log.prior_temp_2 <- scale(log(data$prior_temp_2))
data$z.log.prior_temp_3 <- scale(log(data$prior_temp_3))
data$z.log.prior_temp_4 <- scale(log(data$prior_temp_4))

#Log temp
data$z.log.temp <- scale(log(data$incb_temp))

write.csv(data, row.names = F, "data/data_final/mr_final_analysis.csv")





