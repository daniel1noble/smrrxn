setwd("~/Dropbox/smrrxn/")

rm(list = ls())

library(ggplot2)
library(dplyr)

#data <- read.csv("data/data_final/ldeli_integrate.csv")
#str(data)
#data$date <- as.Date(data$date, "%d/%m/%Y")

#data2 <- data

#correcting bad baselines in CO2 and subtracting CO2 control samples form air samples
#badBLrows <- row.names(data2[which(data2$co2_samp_control > data2$co2_samp_1  | data2$co2_samp_control > data2$co2_samp_2 ),18:20]) 
#length(badBLrows) #64

#data2[rownames(data2) %in% badBLrows,]

#data2[which( data2$co2_samp_control > data2$co2_samp_1  | data2$co2_samp_control > data2$co2_samp_2 ),18] <- NA
 
#data2$co2_samp_1_correct <- ifelse(is.na(data2$co2_samp_control), data2$co2_samp_1,data2$co2_samp_1 - #data2$co2_samp_control)

#data2$co2_samp_2_correct <- ifelse(is.na(data2$co2_samp_control), data2$co2_samp_2,data2$co2_samp_2 - #data2$co2_samp_control)

#data2$frac_co2 <- pmax(data2$co2_samp_1_correct, data2$co2_samp_2_correct)

#write.csv(data2, row.names = F, "data/data_final/corrected_ldeli_integrate.csv") 

#### At this point, I opened the excel file and had to manual input frac_co2 values if there was a NA in either data2$co2_samp_1_correct, data2$co2_samp_2_correct. This is why I read it in R below

finaldata <- read.csv("data/data_final/archive/corrected_ldeli_integrate.csv")
str(finaldata) #2520 obs

#creating series variable
finaldata$series <-paste(finaldata$id , finaldata$samp_period, sep = "_")

#adding population into datafile 

pop.data <- read.csv("data/data_final/male_population.csv")

pop.data <- pop.data[,c(1,4)]
names(pop.data) <- c("id", "site")

newdata <- merge(finaldata, pop.data)
str(newdata) #2541 obs

#obs variable

newdata$obs <- seq(1:nrow(newdata))

#Creating mass total air and total co2 and co2_pmin variables
newdata$lizmass <- newdata$ch_lizmass - newdata$ch_mass
newdata$total_air <- newdata$ch_vol - newdata$lizmass
newdata$total_co2 <- newdata$frac_co2 * newdata$total_air
newdata$co2_pmin <- newdata$total_co2 / newdata$t_diff

#outlier checking

#Potential outliers for co2_pmin:  770, 1006, 1297, 244, 1012, 1483, 1181, 1502, 

#cleveland dot plot for liz mass

ggplot(newdata, aes(lizmass, obs, label = obs)) +
  geom_point(colour = "white") +
  geom_text() + 
  theme_bw()         

#Potential outliers for mass:  359, 377  - typos in ch_mass

newdata[newdata$obs == 359, "ch_mass"] <- 18.195
newdata[newdata$obs == 377, "ch_mass"] <-  17.824 

newdata$lizmass <- newdata$ch_lizmass - newdata$ch_mass

#Potential outliers for mass >2 : 1649 1654 1997 2003 2441 2456 2461 2477 #4 rows of data
newdata[newdata$lizmass > 2,"obs"] 
great2mass.obs <- c(1649, 1654, 1997, 2003, 2441, 2456, 2461, 2477)
newdata[newdata$obs %in% great2mass.obs,]

#obs 1649 1654 not a typo but inconsisent with ID mean mass
ld0128 <- newdata[newdata$id == "ld0128",]
ld0128[,"lizmass"]
summary(subset(ld0128, lizmass < 2, select = lizmass))
#Putting mean of lizard in place of error? or NA?
newdata[newdata$id == "ld0128" & newdata$obs == c(1649, 1654), "lizmass"] #<- rep(1.137,2)

#obs 1997, 2003 ch_mass typo -  need to recalculate liz mass after fixing these ch_mass errors
newdata[newdata$id == "ld0143" & newdata$obs == 1997, "ch_mass"] <- 18.142
newdata[newdata$id == "ld0143" & newdata$obs == 2003, "ch_mass"] <- 18.142

#obs 2441 2456 ch_mass typo - need to recalculate liz mass after fixing these ch_mass errors
newdata[newdata$id == "ld0171" & newdata$obs == 2441, "ch_mass"] <- 18.104
newdata[newdata$id == "ld0171" & newdata$obs == 2456, "ch_mass"] <- 18.104

#obs 2461, 2477 not a typo but inconsisent with ID mean mass *see above, had issues with this lizard's ch_mass too
ld0171  <- newdata[newdata$id == "ld0171",]
ld0171[,"lizmass"]
summary(subset(ld0171, lizmass < 2, select = lizmass))
#Putting mean of lizard in place of error? or NA?
newdata[newdata$id == "ld0171" & newdata$obs == c(2461, 2477), "lizmass"] #<- rep(1.170,2)

#Potential outliers for mass <0.5 : 77  117  146  166 560  563 1306 1327 1457 1496 1601 1623 2134 2176 2287 2288 2395 2403 #9 rows of data
newdata[newdata$lizmass < 0.5,"obs"] 
less0.5mass.obs <- c(77, 117, 146, 166, 560, 563, 1306, 1327, 1457, 1496, 1601, 1623, 2134, 2176, 2287, 2288, 2395, 2403)
newdata[newdata$obs %in% less0.5mass.obs,]

#obs 77, 117 must be a typo, inconsisent with ID mean mass 
ld0006   <- newdata[newdata$id == "ld0006",]
ld0006 [,"lizmass"]
summary(subset(ld0006, lizmass > 0.5, select = lizmass))
#Putting mean of lizard in place of error? or NA?
newdata[newdata$id == "ld0171" & newdata$obs == c(2461, 2477), "lizmass"] #<- rep(1.035,2)

#obs 146, 166 ch_lizmass typo  - need to recalculate liz mass after fixing these ch_lizmass errors
newdata[newdata$id == "ld0011" & newdata$obs == 146, "ch_lizmass"] <- 19.008
newdata[newdata$id == "ld0011" & newdata$obs == 166, "ch_lizmass"] <- 19.008

#obs 560, 563 ch_lizmass typo  - need to recalculate liz mass after fixing these ch_lizmass errors
newdata[newdata$id == "ld0049" & newdata$obs == 560, "ch_lizmass"] <- 19.208
newdata[newdata$id == "ld0049" & newdata$obs == 563, "ch_lizmass"] <- 19.208

#obs 1306, 1327 ch_lizmass typo  - need to recalculate liz mass after fixing these ch_lizmass errors
newdata[newdata$id == "ld0105" & newdata$obs == 1306, "ch_lizmass"] <- 19.451
newdata[newdata$id == "ld0105" & newdata$obs == 1327, "ch_lizmass"] <- 19.451

#obs 1457, 1496 ch_lizmass typo  - need to recalculate liz mass after fixing these ch_lizmass errors
newdata[newdata$id == "ld0115" & newdata$obs == 1457, "ch_lizmass"] <- 19.589
newdata[newdata$id == "ld0115" & newdata$obs == 1496, "ch_lizmass"] <- 19.589

#obs 1601, 1623 ch_lizmass typo  - need to recalculate liz mass after fixing these ch_lizmass errors
newdata[newdata$id == "ld0119" & newdata$obs == 1601, "ch_lizmass"] <- 19.512
newdata[newdata$id == "ld0119" & newdata$obs == 1623, "ch_lizmass"] <- 19.512

#obs 2134, 2176 ch_lizmass typo  - need to recalculate liz mass after fixing these ch_lizmass errors
newdata[newdata$id == "ld0150" & newdata$obs == 2134, "ch_lizmass"] <- 19.033
newdata[newdata$id == "ld0150" & newdata$obs == 2176, "ch_lizmass"] <- 19.033

#obs 2287 2288 ch_lizmass typo  - need to recalculate liz mass after fixing these ch_lizmass errors
newdata[newdata$id == "ld0154" & newdata$obs == 2287, "ch_lizmass"] <- 19.840
newdata[newdata$id == "ld0154" & newdata$obs == 2288, "ch_lizmass"] <- 19.840

#obs 2395 2403 ch_lizmass typo  - need to recalculate liz mass after fixing these ch_lizmass errors
newdata[newdata$id == "ld0170" & newdata$obs == 2395, "ch_lizmass"] <- 19.504
newdata[newdata$id == "ld0170" & newdata$obs == 2403, "ch_lizmass"] <- 19.504

#Potential outliers for  mass > 1.5 and < 2 : 310  362  365  415 1523 1537 #3 rows of data
newdata[newdata$lizmass > 1.5 & newdata$lizmass < 2,"obs"] 
less2great1.5.obs <- c(310, 362, 365, 415, 1523, 1537)
newdata[newdata$obs %in% less2great1.5.obs,]

#obs 310, 362 - not a typo but inconsisent with ID mean mass
ld0020 <- newdata[newdata$id == "ld0020",]
ld0020[,"lizmass"]
summary(subset(ld0020, lizmass < 1.5 , select = lizmass))
#Putting mean of lizard in place of error? or NA?
newdata[newdata$id == "ld0020" & newdata$obs == 310, "lizmass"] #<- 1.094
newdata[newdata$id == "ld0020" & newdata$obs == 362, "lizmass"] #<- 1.094

#obs 365, 415 - not a typo but inconsisent with ID mean mass
ld0027 <- newdata[newdata$id == "ld0027",]
ld0027[,"lizmass"]
summary(subset(ld0027, lizmass < 1.8 , select = lizmass))
#Putting mean of lizard in place of error? or NA?
newdata[newdata$id == "ld0020" & newdata$obs == 365, "lizmass"] #<- 1.348 
newdata[newdata$id == "ld0020" & newdata$obs == 415, "lizmass"] #<- 1.348 

#obs 1523, 1537 - not a typo but inconsisent with ID mean mass
ld0118 <- newdata[newdata$id == "ld0118",]
ld0118[,"lizmass"]
summary(subset(ld0118, lizmass < 1.7 , select = lizmass))
#Putting mean of lizard in place of error? or NA?
newdata[newdata$id == "ld0118" & newdata$obs == 1523, "lizmass"] #<- 1.326 
newdata[newdata$id == "ld0118" & newdata$obs == 1537, "lizmass"] #<- 1.326 

#Potential outliers for  mass > 0.5 and < 0.8 :  83   96  #1 rows of data
newdata[newdata$lizmass > 0.5 & newdata$lizmass < 0.8,"obs"] 
newdata[newdata$obs == c(83,96),]

#obs  83 96 ch_lizmass typo - need to recalculate liz mass after fixing these ch_lizmass errors
newdata[newdata$id == "ld0006" & newdata$obs == 83, "ch_lizmass"] <- 19.537
newdata[newdata$id == "ld0006" & newdata$obs == 96, "lizmass"] <- 19.537

#newdata$lizmass <- newdata$ch_lizmass - newdata$ch_mass

###cleveland dot plot for co2: 
ggplot(newdata,
       aes(x = incb_temp, y = co2_pmin, group = incb_temp, label = obs)) +
  geom_boxplot() +
  geom_text()

ggplot(newdata, aes(co2_pmin, obs, label = obs)) +
  geom_point(colour = "white") +
  geom_text() + 
  theme_bw()

#Potential outliers for co2 > 0.02: 24   25   32  139  244  252  288  295  432  621  660  664  770  773  985 1006 1012 1022 1070 1131 1133 1213 1289 1297 1454 1483 1502 1522 1762 2356 #30 rows of data using obs 432 as cut off
newdata[newdata$obs == 432,"co2_pmin"]
vec <- newdata[newdata$co2_pmin > 0.02 ,"obs"] 
great0.02co2.obs <- vec[!is.na(vec)]
great0.02co2.obs.dat <- newdata[newdata$obs %in% great0.02co2.obs, ]

###cleveland dot plot for incb_temp: #all is fine

ggplot(newdata, aes(incb_temp, obs, label = obs)) +
  geom_point(colour = "white") +
  geom_text() + 
  theme_bw()

###cleveland dot plot for body_temp: 

ggplot(newdata, aes(body_temp, obs, label = obs)) +
  geom_point(colour = "white") +
  geom_text() + 
  theme_bw()

#Outlier in body temp
newdata[newdata$body_temp <10,'obs']

#obs 1394, 1399 TYPO
newdata[newdata$obs == 1394,"body_temp"] <- 24
newdata[newdata$obs == 1399,"body_temp"] <- 24 

newdata[newdata$body_temp >31,'obs']

#obs 1881 1918 #This is not a typo it was what was recorded
newdata[newdata$obs == 1881,] 
newdata[newdata$obs == 1918,] 

newdata[newdata$id == "ld0140","body_temp"]

newdata[newdata$obs == 1881,"body_temp"] 
newdata[newdata$obs == 1918,"body_temp"] 

#Scatterplot of co2~mass with obs labels
#Outliers 

plot(log(newdata$co2_pmin) ~ log(newdata$lizmass))

ggplot(data, aes(y = log(co2_pmin), x = log(lizmass), label = obs)) +
  geom_point(colour = "white") +
  geom_text() + 
  #scale_x_continuous(limits = c(-4, 3)) + 
  #scale_y_continuous(limits = c(-6, 3)) +
  theme_bw()

#Potential outliers: co2~mass obs less than -6.6 on log co2 scale : 113  160  168  242  326  388  401  801 1360 1638 1827 1862 1986 2094 2296 2384
less6.6logco2mass <- newdata[log(newdata$co2_pmin) < -6.6,"obs"]
less6.6logco2mass.obs <- less6.6logco2mass[!is.na(less6.6logco2mass)]

#Potential outliers: co2~mass obs less than -0.8 on log mass scale :  77  117  146  166  560  563 1306 1327 1457 1496 1601 1623 2134 2176 2287 2288 2395 2403
less1logco2mass <- newdata[log(newdata$lizmass) < -0.8,"obs"]
less1logco2mass.obs <- less1logco2mass[!is.na(less1logco2mass)]

#Recalcuation of key variables
newdata$lizmass <- newdata$ch_lizmass - newdata$ch_mass
newdata$total_air <- newdata$ch_vol - newdata$lizmass
newdata$total_co2 <- newdata$frac_co2 * newdata$total_air
newdata$co2_pmin <- newdata$total_co2 / newdata$t_diff

#Write newdata as csv

write.csv(newdata,row.names = F, "data/data_final/mrrxn_nooutlier.csv")

#standardising and transforming variables

data$z.incb_temp <- scale(data$incb_temp)
data$z.log.temp <- scale(log(data$incb_temp))

data$z.body_temp <- scale(data$body_temp)
data$z.log.body_temp <- scale(log(data$body_temp))

data$log.mass <- log(data$lizmass)
data$z.log.mass <- scale(data$log.mass)

data$log.co2pmin <- log(data$co2_pmin)
data$z.log.co2pmin <- scale(data$log.co2pmin)

#Standardising and transforming prior temp variables 

data$z.prior_temp_3 <- scale(data$prior_temp_3)
data$z.log.prior_temp_3 <- scale(log(data$prior_temp_3))

data$z.prior_temp_4 <- scale(data$prior_temp_4)
data$z.log.prior_temp_4 <- scale(log(data$prior_temp_4))

#Boltzman standardisation of temperature predictors
#myTemperatureInKelvin  <-  1:30 + 273.15
#inverseKT              <-  1 / 8.62e-5 * ((1 / mean(myTemperatureInKelvin)) - (1 / myTemperatureInKelvin))
#where 8.62e-5 is the Boltzmann constant (in eV/K). this transformation ensures that your intercept is also independent of temperature (and not only mass), and corresponds to the metabolic rate at mean(myTemperatureInKelvin).

#incb_temp
data$incb_temp_K <- data$incb_temp + 273.15
data$inverseK_incb_temp <- 1 / 8.62e-5 * ((1 / mean(data$incb_temp_K)) - (1 /data$incb_temp_K))

#body_temp
data$body_temp_K <- data$body_temp + 273.15
data$inverseK_body_temp <- 1 / 8.62e-5 * ((1 / mean(data$body_temp_K)) - (1 /data$body_temp_K))

#prior_temp3, 4

