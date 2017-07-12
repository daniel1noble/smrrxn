setwd("~/Dropbox/smrrxn/")

rm(list = ls())

library(ggplot2)
library(lme4)

#histograms

data <- read.csv("data/data_final/mrrxn_final.csv")
str(data)
names(data)

varibs.excl <- c("tc", "date", "ch_num", "key", "ch_mass", "ch_lizmass", "t_flush", "t_samp", "t_diff", "ch_vol", "co2_samp_control", "co2_samp_1", "co2_samp_2", "air_collect_notes", "fms_notes", "co2_warthog_notes", "co2_samp_1_correct", "co2_samp_2_correct", "frac_co2", "total_air", "total_co2")

varibs.need <- c("obs", "samp_period", "id" ,"batch", "series", "incb_num", "incb_temp_id", "defecate", "incb_temp", "z.incb_temp", "z.log.temp", "incb_temp_K", "inverseK_incb_temp", "body_temp", "z.body_temp", "z.log.body_temp", "body_temp_K", "inverseK_body_temp" , "z.prior_temp1", "z.log.prior_temp1", "prior_temp1_K", "inverseK_prior_temp1", "prior_temp2_K", "inverseK_prior_temp2", "z.prior_temp2", "z.log.prior_temp2", "orig_lizmass", "lizmass_nocombout", "log.mass", "z.log.mass","orig_co2_pmin", "co2pm_nocombout", "log.co2pmin", "z.log.co2pmin")

incl.vars <- names(data) %in% varibs.need
data <- data[incl.vars]
View(data)
str(data)

#Normality
#Temp predictors
hist(data$z.incb_temp) ; shapiro.test(data$z.incb_temp)
hist(data$z.log.temp) ; shapiro.test(data$z.log.temp)
hist(data$inverseK_incb_temp) ; shapiro.test(data$inverseK_incb_temp)

hist(data$z.body_temp) ; shapiro.test(data$z.body_temp)
hist(data$z.log.prior_temp1). ; shapiro.test(data$z.log.prior_temp1)
hist(data$inverseK_body_temp) ; shapiro.test(data$inverseK_body_temp)

hist(data$z.prior_temp1) ; shapiro.test(data$z.prior_temp1)
hist(data$z.log.prior_temp1) ; shapiro.test(data$z.log.prior_temp1)
hist(data$inverseK_prior_temp1) ; shapiro.test(data$inverseK_prior_temp1)

hist(data$z.prior_temp2) ; shapiro.test(data$z.prior_temp2)
hist(data$z.log.prior_temp2) ; shapiro.test(data$z.log.prior_temp2)
hist(data$inverseK_prior_temp2) ; shapiro.test(data$inverseK_prior_temp2)

#Mass predictors
hist(data$log.mass) ; shapiro.test(data$log.mass)
hist(data$z.log.mass) ; shapiro.test(data$z.log.mass)

hist(data$log.co2pmin) ; shapiro.test(data$log.co2pmin)
hist(data$z.log.co2pmin) ; shapiro.test(data$z.log.co2pmin)

#batch and incubation number and defecation effects
#on mr
boxplot(data$z.log.co2pmin ~ data$batch)
boxplot(data$z.log.co2pmin ~ data$incb_num)

data$defecate <- ifelse(is.na(data$defecate), 0,1) 
boxplot(data$z.log.co2pmin ~ data$defecate)

#on mass
boxplot(data$z.log.mass ~ data$batch)
boxplot(data$z.log.mass ~ data$incb_num)
boxplot(data$z.log.mass ~ data$defecate)

#plotting to see if there is non-linear relationship
plot(data$z.log.co2pmin  ~ data$inverseK_incb_temp) #Doesn't look like it

plot(data$z.log.co2pmin  ~ data$z.log.mass) #Doesn't look like it

ggplot(data, aes(x = z.log.mass, y = z.log.co2pmin, label = obs)) +
  geom_point(colour = "white") +
  geom_text() + 
  theme_bw() 

#Collinearity
names(data)
varibs <- c("z.log.co2pmin", "inverseK_incb_temp", "z.log.mass", "inverseK_prior_temp1", "inverseK_prior_temp2")
cor(data[,varibs],use = "complete.obs") 

pairs(data[,c(varibs)])

#Use inverse temp or z-transformed
model.1.1 <- lmer(z.log.co2pmin ~ inverseK_incb_temp + z.log.mass + (1+inverseK_incb_temp|id) + (1+inverseK_incb_temp|series), data = data)
summary(model.1.1)
AIC(model.1.1) #5268.247

model.1.2 <- lmer(z.log.co2pmin ~ z.incb_temp + z.log.mass + (1+z.incb_temp|id) + (1+z.incb_temp|series), data = data)
summary(model.1.2)
AIC(model.1.2) #5271.841

model.1.3 <- lmer(z.log.co2pmin ~ z.log.temp + z.log.mass + (1+z.log.temp|id) + (1+z.log.temp|series), data = data)
summary(model.1.3)
AIC(model.1.3) #5266.752  #No AIC diff in z.log.temp and inverseK temp

#Carry over effects of prior temps using inverse K Temp
model.2.1 <- lmer(z.log.co2pmin ~ inverseK_incb_temp + z.log.mass + inverseK_body_temp + (1+inverseK_incb_temp|id) + (1+inverseK_incb_temp|series), data = data)
summary(model.2.1)
AIC(model.2.1) #5205.507

model.2.2 <- lmer(z.log.co2pmin ~ inverseK_incb_temp + z.log.mass + inverseK_prior_temp1 + (1+inverseK_incb_temp|id) + (1+inverseK_incb_temp|series), data = data)
summary(model.2.2)
AIC(model.2.2) #5203.591
plot(model.2.2)
qqnorm(resid(model.2.2))

model.2.3 <- lmer(z.log.co2pmin ~ inverseK_incb_temp + z.log.mass + inverseK_prior_temp2 + (1+inverseK_incb_temp|id) + (1+inverseK_incb_temp|series), data = data)
summary(model.2.3)
AIC(model.2.3) #5145.008 #This is the best final model 
plot(model.2.3)
qqnorm(resid(model.2.3)) 

#Analysis

#R intercept
#print(VarCorr(model.2.1),comp=c("Variance", "Std.Dev."))
mod2.3 <- data.frame(VarCorr(model.2.3, comp = "Variance"))

mod2.3[4,4] / (mod2.3[4,4] + mod2.3[1,4]) #0.5185752

#R slope

mod2.3[5,4] / (mod2.3[5,4] + mod2.3[2,4]) #0.3543051

#R short term

(mod2.3[4,4] +  mod2.3[1,4]) / (mod2.3[4,4] +  mod2.3[1,4] + mod2.3[7,4]) #0.3050962

#R long term

mod2.3[4,4] / (mod2.3[4,4] +  mod2.3[1,4] + mod2.3[7,4]) #0.1582153













