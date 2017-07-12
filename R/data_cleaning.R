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
str(newdata) #2520 obs

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

newdata[newdata$obs == 379, "ch_mass"] <- 18.195
newdata[newdata$obs == 356, "ch_mass"] <-  17.824 

newdata$lizmass <- newdata$ch_lizmass - newdata$ch_mass

#Potential outliers for mass >2 : 1641 1675 1983 1993 2421 2429 2434 2453
newdata[newdata$lizmass > 2,"obs"] 
great2mass.obs <- c(1641, 1675, 1983, 1993, 2421, 2429, 2434, 2453)
newdata[newdata$obs %in% great2mass.obs,]

#obs 1649 1654 not a typo but inconsisent with ID mean mass
ld0128 <- newdata[newdata$id == "ld0128",]
ld0128[,"lizmass"]
summary(subset(ld0128, lizmass < 2, select = lizmass))
#Putting mean of lizard in place of error? or NA?
ld0128[ld0128$obs == 1641, "lizmass"] #<- rep(1.137,2)
ld0128[ld0128$obs == 1675, "lizmass"] #<- rep(1.137,2)

#obs 1983 1993 ch_mass typo -  need to recalculate liz mass after fixing these ch_mass errors
newdata[newdata$id == "ld0143" & newdata$obs == 1983, "ch_mass"] <- 18.142
newdata[newdata$id == "ld0143" & newdata$obs == 1993, "ch_mass"] <- 18.142

#obs 2429 2434 ch_mass typo - need to recalculate liz mass after fixing these ch_mass errors
newdata[newdata$id == "ld0171" & newdata$obs == 2434, "ch_mass"] <- 18.104
newdata[newdata$id == "ld0171" & newdata$obs == 2429, "ch_mass"] <- 18.104

#obs 2421, 2453 not a typo but inconsisent with ID mean mass *see above, had issues with this lizard's ch_mass too
ld0171  <- newdata[newdata$id == "ld0171",]
ld0171[,"lizmass"]
summary(subset(ld0171, lizmass < 2, select = lizmass))
#Putting mean of lizard in place of error? or NA?
ld0171[ld0171$obs == 2421, "lizmass"] #<- rep(1.137,2)
ld0171[ld0171$obs == 2453, "lizmass"] #<- rep(1.137,2)

#Potential outliers for mass <0.5 : 76  112  147  167 560  573 1286 1312 1445 1484 1588 1614 2130 2155 2241 2272 2377 2383
newdata[newdata$lizmass < 0.5,"obs"] 
less0.5mass.obs <- c(76, 112, 147, 167, 560, 573, 1286, 1312, 1445, 1484, 1588, 1614, 2130, 2155, 2241, 2272, 2377, 2383)
newdata[newdata$obs %in% less0.5mass.obs,]

#obs 76, 112 must be a typo, inconsisent with ID mean mass 
ld0006   <- newdata[newdata$id == "ld0006",]
ld0006[,"lizmass"]
summary(subset(ld0006, lizmass > 0.5, select = lizmass))
#Putting mean of lizard in place of error? or NA?
ld0006[ld0006$obs == 76, "lizmass"] #<- rep(1.042,2)
ld0006[ld0006$obs == 112, "lizmass"] #<- rep(1.042,2)

#obs 147  167 ch_lizmass typo  - need to recalculate liz mass after fixing these ch_lizmass errors
newdata[newdata$id == "ld0011" & newdata$obs == 147, "ch_lizmass"] <- 19.008
newdata[newdata$id == "ld0011" & newdata$obs == 167, "ch_lizmass"] <- 19.008

#obs 560, 573 ch_lizmass typo  - need to recalculate liz mass after fixing these ch_lizmass errors
newdata[newdata$id == "ld0049" & newdata$obs == 560, "ch_lizmass"] <- 19.208
newdata[newdata$id == "ld0049" & newdata$obs == 573, "ch_lizmass"] <- 19.208

#obs 1286, 1312 ch_lizmass typo  - need to recalculate liz mass after fixing these ch_lizmass errors
newdata[newdata$id == "ld0105" & newdata$obs == 1286, "ch_lizmass"] <- 19.451
newdata[newdata$id == "ld0105" & newdata$obs == 1312, "ch_lizmass"] <- 19.451

#obs 1445, 1484 ch_lizmass typo  - need to recalculate liz mass after fixing these ch_lizmass errors
newdata[newdata$id == "ld0115" & newdata$obs == 1445, "ch_lizmass"] <- 19.589
newdata[newdata$id == "ld0115" & newdata$obs == 1484, "ch_lizmass"] <- 19.589

#obs 1588, 1614 ch_lizmass typo  - need to recalculate liz mass after fixing these ch_lizmass errors
newdata[newdata$id == "ld0119" & newdata$obs == 1588, "ch_lizmass"] <- 19.512
newdata[newdata$id == "ld0119" & newdata$obs == 1614, "ch_lizmass"] <- 19.512

#obs 2130, 2155 ch_lizmass typo  - need to recalculate liz mass after fixing these ch_lizmass errors
newdata[newdata$id == "ld0150" & newdata$obs == 2130, "ch_lizmass"] <- 19.033
newdata[newdata$id == "ld0150" & newdata$obs == 2155, "ch_lizmass"] <- 19.033

#obs  2241, 2272 ch_lizmass typo  - need to recalculate liz mass after fixing these ch_lizmass errors
newdata[newdata$id == "ld0154" & newdata$obs == 2241, "ch_lizmass"] <- 19.840
newdata[newdata$id == "ld0154" & newdata$obs == 2272, "ch_lizmass"] <- 19.840

#obs 2377, 2383 ch_lizmass typo  - need to recalculate liz mass after fixing these ch_lizmass errors
newdata[newdata$id == "ld0170" & newdata$obs == 2377, "ch_lizmass"] <- 19.504
newdata[newdata$id == "ld0170" & newdata$obs == 2383, "ch_lizmass"] <- 19.504

#Potential outliers for  mass > 1.6 and < 2 : 342  350  414  417 1519 1547 
newdata[newdata$lizmass > 1.6 & newdata$lizmass < 2,"obs"] 
less2great1.5.obs <- c(340, 350, 414, 417, 1519, 1547)
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

#Potential outliers for  mass > 0.5 and < 0.7 :  78   83
newdata[newdata$lizmass > 0.5 & newdata$lizmass < 0.80,"obs"] 
newdata[newdata$obs %in% c(78, 83) ,]

#obs  83 96 ch_lizmass typo - need to recalculate liz mass after fixing these ch_lizmass errors
newdata[newdata$id == "ld0006" & newdata$obs == 78, "ch_lizmass"] <- 19.537
newdata[newdata$id == "ld0006" & newdata$obs == 83, "ch_lizmass"] <- 19.537

#newdata$lizmass <- newdata$ch_lizmass - newdata$ch_mass
#Combining all true mass outliers together

truemassout.obs <- c(1641, 1675, 2421, 2453, 76, 112, 310, 362, 365, 415, 1523, 1537)
truemassout <- newdata[newdata$obs %in% truemassout.obs,]

newdata$mass_outlier <- ifelse(newdata$obs %in% truemassout.obs, 1,  0) 

###cleveland dot plot for co2: 
ggplot(newdata,
       aes(x = incb_temp, y = co2_pmin, group = incb_temp, label = obs)) +
  geom_boxplot() +
  geom_text()

ggplot(newdata, aes(co2_pmin, obs, label = obs)) +
  geom_point(colour = "white") +
  geom_text() + 
  theme_bw()

#Potential outliers for co2 > 0.02328337 (obs 2304 as cutoff): 249  766 1015 1298
newdata[newdata$obs == 2304, "co2_pmin"]
vec <- newdata[newdata$co2_pmin > 0.02 ,"obs"] 
great0.02co2.obs <- vec[!is.na(vec)]
co2outs <- newdata[newdata$obs %in% great0.02co2.obs, ]
write.csv(co2outs , row.names = F, "data/co2outliers.csv")

#Treat the ones with notes in aircollect, fms and warthog as outliers as they were mechanic errors 

co2trueouts.obs <- c(979, 1053, 1015, 262, 1273, 274, 633, 1122, 1443, 1298 ,249 ,766, 1230, 645, 1471, 1490, 994, 2304, 1530, 14)
co2trueouts <- newdata[newdata$obs %in% co2trueouts.obs, ]
co2trueouts[order(co2trueouts$obs),]

newdata$co2_outlier <- ifelse(newdata$obs %in% co2trueouts.obs, 1,  0) 
 
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

#obs 1387 1390 TYPO
newdata[newdata$obs == 1387,"body_temp"] <- 24
newdata[newdata$obs == 1390,"body_temp"] <- 24 

newdata[newdata$body_temp >31,'obs']

#obs 1881 1918 #This is not a typo it was what was recorded - decide to leave this
newdata[newdata$obs == 1864,] 
newdata[newdata$obs == 1892,] 

newdata[newdata$id == "ld0140","body_temp"]

#Scatterplot of co2~mass with obs labels
#Outliers 

plot(log(newdata$co2_pmin) ~ log(newdata$lizmass))

#Potential outliers: co2~mass obs less than -7.5 on log co2 scale : 1348, 1836
less7.5logco2mass <- newdata[log(newdata$co2_pmin) < -7.5,"obs"]
less7.5logco2mass.obs <- less7.5logco2mass[!is.na(less7.5logco2mass)]

#Potential outliers: co2~mass obs less than -0.5 on log mass scale :  76 112
less0.5logco2mass <- newdata[log(newdata$lizmass) < -0.5,'obs']
less0.5logco2mass.obs <- less0.5logco2mass[!is.na(less0.5logco2mass)]

#Potential outliers: co2~mass obs greater than 0.5 on log mass scale :  414  417 1519 1547 1641 1675 2421 2453
great0.5logco2mass <- newdata[log(newdata$lizmass) > 0.5,'obs']
great0.5logco2mass.obs <- great0.5logco2mass[!is.na(great0.5logco2mass)]

co2massout <- c(less7.5logco2mass.obs, less0.5logco2mass.obs, great0.5logco2mass)
co2massout.obs <- co2massout[!is.na(co2massout)]

newdata$co2mass_outlier <- ifelse(newdata$obs %in% co2massout.obs, 1,  0)

#Boltzman standardisation of temperature predictors
#myTemperatureInKelvin  <-  1:30 + 273.15
#inverseKT              <-  1 / 8.62e-5 * ((1 / mean(myTemperatureInKelvin)) - (1 / myTemperatureInKelvin))
#where 8.62e-5 is the Boltzmann constant (in eV/K). this transformation ensures that your intercept is also independent of temperature (and not only mass), and corresponds to the metabolic rate at mean(myTemperatureInKelvin).

#incb_temp
newdata$incb_temp_K <- newdata$incb_temp + 273.15
newdata$inverseK_incb_temp <- 1 / 8.62e-5 * ((1 / mean(newdata$incb_temp_K)) - (1 /newdata$incb_temp_K))

#body_temp
newdata$body_temp_K <- newdata$body_temp + 273.15
newdata$inverseK_body_temp <- 1 / 8.62e-5 * ((1 / mean(newdata$body_temp_K, na.rm =T)) - (1 /newdata$body_temp_K))

#Write newdata as csv

write.csv(newdata,row.names = F, "data/data_final/mrrxn_orig.csv")

#Read in data again after putting NA in mass and co2 if allcomb_outlier indicator = 1 and creating prior_temp1 and prior_temp2
newdata <- read.csv("data/data_final/mrrxn_orig.csv")
str(newdata) #2520 obs

#Recalcuation of key variables
newdata$lizmass <- newdata$ch_lizmass - newdata$ch_mass
newdata$total_air <- newdata$ch_vol - newdata$lizmass
newdata$total_co2 <- newdata$frac_co2 * newdata$total_air
newdata$co2_pmin <- newdata$total_co2 / newdata$t_diff

#standardising and transforming variables

newdata$z.incb_temp <- scale(newdata$incb_temp)
newdata$z.log.temp <- scale(log(newdata$incb_temp))

newdata$z.body_temp <- scale(newdata$body_temp)
newdata$z.log.body_temp <- scale(log(newdata$body_temp))

newdata$log.mass <- log(newdata$lizmass_nocombout)
newdata$z.log.mass <- scale(newdata$log.mass)

newdata$log.co2pmin <- log(newdata$co2pm_nocombout)
newdata$z.log.co2pmin <- scale(newdata$log.co2pmin)

#Standardising and temp correctiong prior temp variables

#incb_temp
newdata$incb_temp_K <- newdata$incb_temp + 273.15
newdata$inverseK_incb_temp <- 1 / 8.62e-5 * ((1 / mean(newdata$incb_temp_K)) - (1 /newdata$incb_temp_K))

#body_temp
newdata$body_temp_K <- newdata$body_temp + 273.15
newdata$inverseK_body_temp <- 1 / 8.62e-5 * ((1 / mean(newdata$body_temp_K, na.rm =T)) - (1 /newdata$body_temp_K))

newdata$z.prior_temp1 <- scale(newdata$prior_temp1)
newdata$z.log.prior_temp1 <- scale(log(newdata$prior_temp1))

newdata$z.prior_temp2 <- scale(newdata$prior_temp2)
newdata$z.log.prior_temp2 <- scale(log(newdata$prior_temp2))

#prior_temp 1
newdata$prior_temp1_K <- newdata$prior_temp1 + 273.15
newdata$inverseK_prior_temp1 <- 1 / 8.62e-5 * ((1 / mean(newdata$prior_temp1_K, na.rm =T)) - (1 /newdata$prior_temp1_K))

#prior_temp 2
newdata$prior_temp2_K <- newdata$prior_temp2 + 273.15
newdata$inverseK_prior_temp2 <- 1 / 8.62e-5 * ((1 / mean(newdata$prior_temp2_K, na.rm =T)) - (1 /newdata$prior_temp2_K))

View(newdata)
write.csv(newdata, row.names = F, "data/data_final/mrrxn_final.csv")
