setwd("~/Dropbox/smrrxn/")

rm(list = ls())

library(ggplot2)
library(lme4)

#histograms

data <- read.csv("data/data_final/mr_final_analysis.csv")

#Normality

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

#plotting to see if there is non-linear relationship
boxplot(data$co2_pmin ~ data$incb_temp) #Doesn't look like it

plot(log(data$co2_pmin) ~ log(data$lizmass)) #Doesn't look like it

#Carryover effects of temperature order
model.1.1 <- lmer(z.log.co2pmin ~ z.incb_temp + z.log.mass + z.body_temp + (1+z.incb_temp|id) + (1+z.incb_temp|series), data = data)
summary(model.1.1)
AIC(model.1.1) #5341.076

model.1.2 <- lmer(z.log.co2pmin ~ z.incb_temp + z.log.mass + z.prior_temp_3 + (1+z.incb_temp|id) + (1+z.incb_temp|series), data = data)
summary(model.1.2)
AIC(model.1.2) #5345.327

plot(data$co2_pmin ~ data$z.prior_temp_3)
plot(data$co2_pmin ~ data$z.log.prior_temp_3)

model.1.3 <- lmer(z.log.co2pmin ~ z.incb_temp + z.log.mass + z.log.prior_temp_4 + (1+z.incb_temp|id) + (1+z.incb_temp|series), data = data)
summary(model.1.3)
AIC(model.1.3) #5278.598

plot(data$co2_pmin ~ data$z.prior_temp_4)
plot(data$co2_pmin ~ data$z.log.prior_temp_4)

#Collinearity
names(data)
varibs <- c("z.log.co2pmin", "incb_temp", "z.log.mass", "z.prior_temp_3", "z.prior_temp_4")
cor(data[,varibs],use = "complete.obs") 

pairs(data[,c(varibs)], log = "xy")

#Analysis
#To use log temp or not
model.2 <- lmer(z.log.co2pmin ~ z.incb_temp + z.log.mass + z.prior_temp_4 + (1+z.incb_temp|id) + (1+z.incb_temp|series), data = data)
summary(model.2)
AIC(model.2) #5278.598
plot(model.2)
qqnorm(resid(model.2))

model.2.2 <- lmer(z.log.co2pmin ~ z.log.temp + z.log.mass + z.prior_temp_4 + (1+z.log.temp|id) + (1+z.log.temp|series), data = data)
summary(model.2.2)
AIC(model.2.2) #5267.552 #Only a little better using log temp
plot(model.2.2)
qqnorm(resid(model.2.2))

#R intercept
#print(VarCorr(model.2.1),comp=c("Variance", "Std.Dev."))
mod2.2 <- data.frame(VarCorr(model.2.2, comp = "Variance"))

mod2.2[4,4] / (mod2.2[4,4] + mod2.2[1,4]) #0.494933

#R slope

mod2.2[5,4] / (mod2.2[5,4] + mod2.2[2,4]) #0.2683569

#R short term

(mod2.2[4,4] +  mod2.2[1,4]) / (mod2.2[4,4] +  mod2.2[1,4] + mod2.2[7,4])

#R long term

mod2.2[4,4] / (mod2.2[4,4] +  mod2.2[1,4] + mod2.2[7,4])













