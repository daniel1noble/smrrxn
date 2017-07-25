setwd("~/Dropbox/smrrxn/")

rm(list = ls())

library(ggplot2)
library(lme4)

data <- read.csv("data/data_final/mrrxn_final_v2.csv")
str(data)
names(data)

#varibs.excl <- c("tc", "date", "ch_num", "key", "ch_mass", "ch_lizmass", "t_flush", "t_samp", "t_diff", "ch_vol", "co2_samp_control", "co2_samp_1", "co2_samp_2", "air_collect_notes", "fms_notes", "co2_warthog_notes", "co2_samp_1_correct", "co2_samp_2_correct", "frac_co2", "total_air", "total_co2")

varibs.need <- c("obs", "samp_period", "id" ,"batch", "series", "incb_num", "incb_temp_id", "defecate", "incb_temp", "z.incb_temp", "z.log.temp", "incb_temp_K", "inverseK_incb_temp", "body_temp", "z.body_temp", "z.log.body_temp", "body_temp_K", "inverseK_body_temp" , "z.prior_temp1", "z.log.prior_temp1", "prior_temp1_K", "inverseK_prior_temp1", "prior_temp2_K", "inverseK_prior_temp2", "z.prior_temp2", "z.log.prior_temp2", "orig_lizmass", "lizmass_nocombout", "log.mass", "z.log.mass","orig_co2_pmin", "co2pm_nocombout", "log.co2pmin", "z.log.co2pmin")

incl.vars <- names(data) %in% varibs.need
data <- data[incl.vars]
#View(data)
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

#plotting to see if there is non-linear relationship with temperature
plot(data$z.log.co2pmin  ~ data$inverseK_incb_temp) #Doesn't look like it

plot(data$z.log.co2pmin  ~ data$z.log.mass, ylim=c(3,-3)) #Doesn't look like it

scatter.smooth(data$z.log.co2pmin ~ data$inverseK_incb_temp) 
scatter.smooth(data$z.log.co2pmin ~ data$inverseK_incb_temp) 

ggplot(data, aes(y = z.log.co2pmin, x = z.log.mass, obs, label = obs)) +
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
AIC(model.2.3) #5106.399 #This is the best final model 
plot(model.2.3)
qqnorm(resid(model.2.3)) 

#Check if my slope and intercept for mass is similar to Uyedas
model.2.4 <- lmer(log.co2pmin ~ inverseK_incb_temp + log.mass + inverseK_prior_temp2 + (1+inverseK_incb_temp|id) + (1+inverseK_incb_temp|series), data = data)
summary(model.2.4)

#Checking when temperature and body mass is adding variation to within or between individual variation? 
model.2.5 <- lmer(log.co2pmin ~ (1|id) + (1|series), data = data)
summary(model.2.5)
mod2.5 <- data.frame(VarCorr(model.2.5, comp = "Variance"))
(mod2.5[2,4] + mod2.5[1,4]) / (mod2.5[2,4] + mod2.5[1,4] +  mod2.5[3,4])

model.2.6 <- lmer(log.co2pmin ~ inverseK_incb_temp + (1|id) + (1|series), data = data)
summary(model.2.6)
mod2.6 <- data.frame(VarCorr(model.2.6, comp = "Variance"))
(mod2.6[2,4] + mod2.6[1,4]) / (mod2.6[2,4] + mod2.6[1,4] +  mod2.6[3,4])

model.2.7 <- lmer(log.co2pmin ~ inverseK_incb_temp + z.log.mass + (1|id) + (1|series), data = data)
summary(model.2.7)
mod2.7 <- data.frame(VarCorr(model.2.7, comp = "Variance"))
(mod2.7[2,4] + mod2.7[1,4]) / (mod2.7[2,4] + mod2.7[1,4] +  mod2.7[3,4])

#Analysis Repeatability of metabolic rate

#R intercept
#print(VarCorr(model.2.1),comp=c("Variance", "Std.Dev."))
mod2.3 <- data.frame(VarCorr(model.2.3, comp = "Variance"))

mod2.3[4,4] / (mod2.3[4,4] + mod2.3[1,4]) 

#R slope

mod2.3[5,4] / (mod2.3[5,4] + mod2.3[2,4]) 

#R short term

(mod2.3[4,4] +  mod2.3[1,4]) / (mod2.3[4,4] +  mod2.3[1,4] + mod2.3[7,4]) 
#R long term

mod2.3[4,4] / (mod2.3[4,4] +  mod2.3[1,4] + mod2.3[7,4])

#Repeatabilty of metabolic rate at each temp
predictors <- c("z.log.co2pmin", "z.log.mass", "inverseK_prior_temp2", "id", "series")

unique(data$incb_temp)

twtwo_dat <- data[data$incb_temp == 22,]
model.3.1 <- lmer(z.log.co2pmin ~ z.log.mass + inverseK_prior_temp2 + samp_period + (1|id), data = twtwo_dat)
summary(model.3.1)
plot(model.3.1)
qqnorm(resid(model.3.1)) 
mod3.1 <- data.frame(VarCorr(model.3.1, comp = "Variance"))
mod3.1[1,4] / (mod3.1[2,4] + mod3.1[1,4]) 

twfr_dat <- data[data$incb_temp == 24,]
model.3.2 <- lmer(z.log.co2pmin ~ z.log.mass + inverseK_prior_temp2 + samp_period + (1|id), data = twfr_dat)
summary(model.3.2)
plot(model.3.2)
qqnorm(resid(model.3.2))
mod3.2 <- data.frame(VarCorr(model.3.2, comp = "Variance"))
mod3.2[1,4] / (mod3.2[2,4] + mod3.2[1,4]) 

twsx_dat <- data[data$incb_temp == 26,]
model.3.3 <- lmer(z.log.co2pmin ~ z.log.mass + inverseK_prior_temp2 + samp_period + (1|id), data = twsx_dat)
summary(model.3.3)
plot(model.3.3)
qqnorm(resid(model.3.3))
mod3.3 <- data.frame(VarCorr(model.3.3, comp = "Variance"))
mod3.3[1,4] / (mod3.3[2,4] + mod3.3[1,4]) 

twet_dat <- data[data$incb_temp == 28,]
model.3.4 <- lmer(z.log.co2pmin ~ z.log.mass + inverseK_prior_temp2 + samp_period + (1|id), data = twet_dat)
summary(model.3.4)
plot(model.3.4)
qqnorm(resid(model.3.4))
mod3.4 <- data.frame(VarCorr(model.3.4, comp = "Variance"))
mod3.4[1,4] / (mod3.4[2,4] + mod3.4[1,4]) 

thrt_dat <- data[data$incb_temp == 30,]
model.3.5 <- lmer(z.log.co2pmin ~ z.log.mass + inverseK_prior_temp2 + samp_period + (1|id), data = thrt_dat)
summary(model.3.5)
plot(model.3.5)
qqnorm(resid(model.3.5))
mod3.5 <- data.frame(VarCorr(model.3.5, comp = "Variance"))
mod3.5[1,4] / (mod3.5[2,4] + mod3.5[1,4]) 

thrttw_dat <- data[data$incb_temp == 32,]
model.3.6 <- lmer(z.log.co2pmin ~ z.log.mass + inverseK_prior_temp2 + samp_period + (1|id), data = thrttw_dat)
summary(model.3.6)
plot(model.3.6)
qqnorm(resid(model.3.6))
mod3.6 <- data.frame(VarCorr(model.3.6, comp = "Variance"))
mod3.6[1,4] / (mod3.6[2,4] + mod3.6[1,4]) 

rpt.plot <- data.frame(temp = sort(unique(data$incb_temp)),
              rpt = c( mod3.1[1,4] / (mod3.1[2,4] + mod3.1[1,4]),
                       mod3.2[1,4] / (mod3.2[2,4] + mod3.2[1,4]),
                       mod3.3[1,4] / (mod3.3[2,4] + mod3.3[1,4]),
                       mod3.4[1,4] / (mod3.4[2,4] + mod3.4[1,4]),
                       mod3.5[1,4] / (mod3.5[2,4] + mod3.5[1,4]),
                       mod3.6[1,4] / (mod3.6[2,4] + mod3.6[1,4])))

ggplot(rpt.plot, aes(y = rpt,
                     x = temp)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(name = "Temperature",
                     limits= c(22, 32),
                     breaks= seq(22,32,2)) + 
  scale_y_continuous(name = "Adjusted repeatability",
                     limits = c(0, 0.25)) + 
  theme_bw()   

#Trying to plot these reaction norms out by groups of 14 lizards

ids_order <- sort(unique(data$id))
g1 <- ids_order[1:14]
g2 <- ids_order[15:28]
g3 <- ids_order[29:42]

predictors <- c("z.log.co2pmin", "inverseK_incb_temp", "z.log.mass", "inverseK_prior_temp2" , "id" ,"series")

plot.data <- data[complete.cases(data[ , predictors ]),]
plot.data$z.log.co2pm_pred <- predict(model.2.3, type = "response")
plot.data$log.co2pm_pred <- (plot.data$z.log.co2pm_pred * sd(plot.data$log.co2pmin)) + mean(plot.data$log.co2pmin)
plot.data$co2pm_pred <- exp(plot.data$log.co2pm_pred)

g1.dat <- plot.data[plot.data$id %in% g1,]
g2.dat <- plot.data[plot.data$id %in% g2,]
g3.dat <- plot.data[plot.data$id %in% g3,]

#g1
g1.dat$plot_samp_p <- as.factor(g1.dat$samp_period)
levels(g1.dat$plot_samp_p) <- paste0("Sampling series ", levels(g1.dat$plot_samp_p))

ggplot(g1.dat, aes(y = z.log.co2pm_pred, 
                   x = inverseK_incb_temp,
                   group = id)) +
  geom_point(shape = 1, fill = "white", size = 1) +
  geom_line(aes(colour = id)) + 
  facet_wrap( ~ plot_samp_p, nrow = 2) +
  labs(x = "Temperature (1/K)", y = expression(Metabolic~rate~(CO[2]~min^{-1}))) +
  theme_bw() + 
  theme(legend.position = "none",
      #panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"))

ggplot(g1.dat, aes(y = z.log.co2pm_pred, 
                   x = inverseK_incb_temp,
                   group = id, 
                   colour = id)) +
  geom_point() +
  geom_line() + 
  facet_wrap( ~ samp_period, nrow = 2) +
  theme_bw()    

#g2
ggplot(g2.dat, aes(y = z.log.co2pm_pred, 
                   x = inverseK_incb_temp,
                   group = id, 
                   colour = id)) +
  geom_point() +
  geom_line() + 
  facet_wrap( ~ samp_period, nrow = 2) +
  theme_bw() 

#g3
ggplot(g3.dat, aes(y = z.log.co2pm_pred, 
                   x = inverseK_incb_temp,
                   group = id, 
                   colour = id)) +
  geom_point() +
  geom_line() + 
  facet_wrap( ~ samp_period, nrow = 2) +
  theme_bw() 

#Plot all reaction norms for each series, for each individual 
ggplot(g1.dat, aes(y = z.log.co2pm_pred, 
                   x = inverseK_incb_temp,
                   group = as.factor(samp_period), 
                   colour = as.factor(samp_period))) +
  geom_point() +
  geom_line() + 
  facet_wrap( ~ id, nrow = 2) +
  theme_bw() +
  scale_colour_discrete(name = "Sampling period") 

ggplot(g2.dat, aes(y = z.log.co2pm_pred, 
                   x = inverseK_incb_temp,
                   group = as.factor(samp_period), 
                   colour = as.factor(samp_period))) +
  geom_point() +
  geom_line() + 
  facet_wrap( ~ id, nrow = 2) +
  theme_bw() +
  scale_colour_discrete(name = "Sampling period") 

ggplot(g3.dat, aes(y = z.log.co2pm_pred, 
                   x = inverseK_incb_temp,
                   group = as.factor(samp_period), 
                   colour = as.factor(samp_period))) +
  geom_point() +
  geom_line() + 
  facet_wrap( ~ id, nrow = 2) +
  theme_bw() +
  scale_colour_discrete(name = "Sampling period") 

#Try and plot intercept and slope for each lizard

intslope.plot <- coef(model.2.3)$series
nrow(intslope.plot)

intslope.plot$series <- row.names(intslope.plot)
unlist_series <- unlist(strsplit(intslope.plot$series, "_"))

write.csv(intslope.plot, "test.csv")
intslope.dat <- read.csv("data/intslope.csv")

ggplot(intslope.dat, aes(y = Intercept, 
                         x = inverseK_incb_temp, 
                         colour = id)) +
  geom_point() +
  scale_x_continuous(name = "Slope for temperature") + 
  facet_wrap( ~ samp_period, nrow = 2) +
  theme_bw() +
  theme(legend.position = "none")

#Calculating R^2c and R^2m and PCV
#Full model
summary(model.2.3)

#Intercept model
model.2.3.0 <- lmer(z.log.co2pmin ~ 1 + (1+inverseK_incb_temp|id) + (1+inverseK_incb_temp|series), data = data)
summary(model.2.3.0)

# Extraction of fitted value for the alternative model
# fixef() extracts coefficents for fixed effects
# mF@pp$X returns fixed effect design matrix

Fixed <- fixef(model.2.3)[2] * model.2.3@pp$X[, 2] + fixef(model.2.3)[3] * model.2.3@pp$X[, 3] + fixef(model.2.3)[4] * model.2.3@pp$X[, 4]

# Calculation of the variance in fitted values
VarF <- var(Fixed)

# An alternative way for getting the same result
VarF <- var(as.vector(fixef(model.2.3) %*% t(model.2.3@pp$X)))

# R2GLMM(m) - marginal R2GLMM i.e. just fixed effects 
# Equ. 26, 29 and 30
# VarCorr() extracts variance components
# attr(VarCorr(lmer.model),'sc')^2 extracts the residual variance
# marginal R2GLMM i.e. just fixed effects for RANDOM INTERCEPTS ONLY
VarF/(VarF + VarCorr(model.2.3)$id[1,1]  + VarCorr(model.2.3)$series[1,1] + attr(VarCorr(model.2.3), "sc")^2)

# marginal R2GLMM i.e. just fixed effects for RANDOM SLOPES
VarF/(VarF  + VarCorr(model.2.3)$id[2,2]  + VarCorr(model.2.3)$series[2,2] + attr(VarCorr(model.2.3), "sc")^2)

# marginal R2GLMM i.e. just fixed effects for RANDOM INTERCEPTS AND SLOPES
VarF/(VarF + VarCorr(model.2.3)$id[1,1] + VarCorr(model.2.3)$id[2,2] + VarCorr(model.2.3)$series[1,1] + VarCorr(model.2.3)$series[2,2]+ attr(VarCorr(model.2.3), "sc")^2)


# R2GLMM(c) - conditional R2GLMM for for RANDOM INTERCEPTS ONLY
(VarF + VarCorr(model.2.3)$id[1,1]  + VarCorr(model.2.3)$series[1,1])/(VarF + VarCorr(model.2.3)$id[1,1]  + VarCorr(model.2.3)$series[1,1] + attr(VarCorr(model.2.3), "sc")^2)

# R2GLMM(c) - conditional R2GLMM for for RANDOM SLOPES ONLY
(VarF + VarCorr(model.2.3)$id[2,2]  + VarCorr(model.2.3)$series[2,2])/(VarF + VarCorr(model.2.3)$id[2,2]  + VarCorr(model.2.3)$series[2,2] + attr(VarCorr(model.2.3), "sc")^2)

# R2GLMM(c) - conditional R2GLMM for for RANDOM INTERCEPTS AND SLOPES
(VarF + VarCorr(model.2.3)$id[1,1]  + VarCorr(model.2.3)$series[1,1] + VarCorr(model.2.3)$id[2,2]  + VarCorr(model.2.3)$series[2,2])/(VarF + VarCorr(model.2.3)$id[1,1]  + VarCorr(model.2.3)$series[1,1] + VarCorr(model.2.3)$id[2,2]  + VarCorr(model.2.3)$series[2,2] + attr(VarCorr(model.2.3), "sc")^2)


# PCVid0 - proportional change in variance at ID intercept level
(1 - (VarCorr(model.2.3)$id[1,1] / VarCorr(model.2.3.0)$id[1,1])) * 100

# PCVid1 - proportional change in variance at ID slope level
(1 - (VarCorr(model.2.3)$id[2,2] / VarCorr(model.2.3.0)$id[2,2])) * 100

# PCVcontainer - proportional change in variance at series intercept level
(1 - (VarCorr(model.2.3)$series[1,1] / VarCorr(model.2.3.0)$series[1,1])) * 100

# PCVcontainer - proportional change in variance at series slope level
(1 - (VarCorr(model.2.3)$series[2,2] / VarCorr(model.2.3.0)$series[2,2])) * 100

# PCVresid - proportional change in variance at units level
# Equ. 33
(1 - (attr(VarCorr(model.2.3), "sc")^2) / (attr(VarCorr(model.2.3.0), "sc")^2)) * 100
