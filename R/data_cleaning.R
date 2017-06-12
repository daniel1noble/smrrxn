setwd("~/Dropbox/smrrxn/")

data <- read.csv("data/data_final/ldeli_integrate.csv")
str(data)
data$date <- as.Date(data$date, "%d/%m/%Y")

data2 <- data

#correcting bad baselines in CO2 and subtracting CO2 control samples form air samples
badBLrows <- row.names(data2[which(data2$co2_samp_control > data2$co2_samp_1  | data2$co2_samp_control > data2$co2_samp_2 ),18:20]) 
length(badBLrows) #64

data2[rownames(data2) %in% badBLrows,]

data2[which( data2$co2_samp_control > data2$co2_samp_1  | data2$co2_samp_control > data2$co2_samp_2 ),18] <- NA
 
data2$co2_samp_1_correct <- ifelse(is.na(data2$co2_samp_control), data2$co2_samp_1,data2$co2_samp_1 - data2$co2_samp_control)

data2$co2_samp_2_correct <- ifelse(is.na(data2$co2_samp_control), data2$co2_samp_2,data2$co2_samp_2 - data2$co2_samp_control)

data2$frac_co2 <- pmax(data2$co2_samp_1_correct, data2$co2_samp_2_correct)

write.csv(data2, row.names = F, "data/data_final/corrected_ldeli_integrate.csv") 

finaldata <- read.csv("data/data_final/corrected_ldeli_integrate.csv")
str(finaldata) #2541 obs

finaldata$lizmass <- finaldata$ch_lizmass - finaldata$ch_mass
finaldata$total_air <- finaldata$ch_vol - finaldata$lizmass
finaldata$total_co2 <- finaldata$frac_co2 * finaldata$total_air
finaldata$co2_pmin <- finaldata$total_co2 / finaldata$t_diff

#creating series variable
finaldata$series <-paste(finaldata$id , finaldata$samp_period, sep = "_")

#adding population into datafile 

pop.data <- read.csv("data/data_final/male_population.csv")

pop.data <- pop.data[,c(1,4)]
names(pop.data) <- c("id", "site")

newdata <- merge(finaldata, pop.data)
str(newdata)

#obs variable

newdata$obs <- seq(1:nrow(newdata))

#data exploring
#outlier checking

#ggplot(data, aes(as.factor(incb_temp), co2_mlpmin)) +
# geom_boxplot()

boxplot(data$co2_mlpmin ~ data$incb_temp)

#cleveland dot plot

#not sure why base plotting is not working
#plot(data$obs_num ~ data$co2_mlpmin)
#text(data$obs_num, data$co2_mlpmin, labels = data$obs_num,cex = 1)


ggplot(newdata, aes(lizmass, obs, label = obs)) +
  geom_point(colour = "white") +
  geom_text(check_overlap = T) + 
  theme_bw()         

#Potential outliers for mass:  359, 377 

newdata[newdata$obs == 359,12] <- 18.195
newdata[newdata$obs == 377,12]<-  17.824 

newdata$lizmass <- newdata$ch_lizmass - newdata$ch_mass
newdata$total_air <- newdata$ch_vol - newdata$lizmass
newdata$total_co2 <- newdata$frac_co2 * newdata$total_air
newdata$co2_pmin <- newdata$total_co2 / newdata$t_diff

#Potential outliers for mass:  2395, 2287, 2134, 1601, 1457, 560, 146

newdata[newdata$obs == 2395,13] <- 19.504
newdata[newdata$obs == 2287,13] <- 19.840
newdata[newdata$obs == 2134,13] <- 19.033
newdata[newdata$obs == 1601,13] <- 19.512
newdata[newdata$obs == 1457,13] <- 19.589
newdata[newdata$obs == 560,13] <- 19.208
newdata[newdata$obs == 146,13] <- 19.008
newdata[newdata$obs == 77,13] <- newdata[newdata$obs == 77,12]  + mean(newdata[newdata$id == 'ld0006',13]) - mean(newdata[newdata$id == 'ld0006',12]) #Writing error, data in datasheet not correct either.

#Potential outliers for mass:  2403, 2288, 2176, 1623, 1496, 1306, 560, 166, 117

newdata[newdata$obs == 2403,13] <- 19.504
newdata[newdata$obs == 2288,13] <- 19.840
newdata[newdata$obs == 2176,13] <- 19.033
newdata[newdata$obs == 1623,13] <- 19.512
newdata[newdata$obs == 1496,13] <- 19.589
newdata[newdata$obs == 1306,13] <- 19.451
newdata[newdata$obs == 560,13] <- 19.208
newdata[newdata$obs == 166,13] <- 19.008
newdata[newdata$obs == 117,13] <- newdata[newdata$obs == 77,13]

#Potential outliers for mass: 1327, 563
newdata[newdata$obs == 1327,13] <- 19.451
newdata[newdata$obs == 563,13] <- 19.208

#Potential outliers for mass: 2441, 1997, 1649, 365, 1523, 310
newdata[newdata$obs == 2441,12] <- 18.104
newdata[newdata$obs == 1997,12] <- 18.142 
newdata[newdata$obs == 1649,13] <- 17.594 + 1.226 #Writing error, used day prior mass
newdata[newdata$obs == 365,13] <- 17.233 + 1.159 #Writing error, used day prior mass
newdata[newdata$obs == 1523,13] <- 17.422 + 1.420 #Writing error, used day prior mass
newdata[newdata$obs == 310,13] <- 17.199 + 1.191
newdata[newdata$obs == 362,13] <- 17.199 + 1.191

#Potential outliers for mass: 2456, 2003, 1654, 1537, 415

newdata[newdata$obs == 2456,12] <- 18.104
newdata[newdata$obs == 2003,12] <- 18.142 
newdata[newdata$obs == 1654,13] <- 17.594 + 1.226 #Writing error, used day prior mass 
newdata[newdata$obs == 1537,13] <- 17.422 + 1.420 #Writing error, used day prior mass
newdata[newdata$obs == 415,13] <- 17.233 + 1.159 #Writing error, used day prior mass 
newdata[newdata$obs == 2461,13] <-  17.849 + 1.214
newdata[newdata$obs == 2477,13] <-  17.849 + 1.214

#Potential outliers for mass: 2301, 83

newdata[newdata$obs == 2301,] #Not an outlier 
newdata[newdata$obs == 83,13] <- 19.537
newdata[newdata$obs == 96,13] <- 19.537

#Potential outliers for co2: 1297, 770, 244

ggplot(newdata, aes(co2_pmin, obs, label = obs)) +
  geom_point(colour = "white") +
  geom_text(check_overlap = T) + 
  theme_bw()

newdata[newdata$obs == 1297,31] <- NA #BadCO2S

newdata[newdata$obs == 770, 31] <- NA  #BigPT1S2? 
newdata[newdata$id == 'ld0069' & newdata$incb_temp == "32",]

newdata[newdata$obs == 244,31] <- NA  #BadT2BL  
newdata[newdata$id == 'ld0018' & newdata$incb_temp == "30",]

newdata[newdata$obs == 1006,] #Could be an outlier
newdata[newdata$id == 'ld0088' & newdata$incb_temp == "32",]

#Write that csv

write.csv(newdata, row.names = F, "data/data_final/smr_rawfinal.csv")


