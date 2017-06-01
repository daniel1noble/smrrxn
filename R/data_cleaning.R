setwd("~/Dropbox/smrrxn/")

data <- read.csv("data/data_final/ldeli_integrate.csv")
str(data)
data$date <- as.Date(data$date, "%d/%m/%Y")

data2 <- data

badBLrows <- row.names(data2[which(data2$co2_samp_control > data2$co2_samp_1  | data2$co2_samp_control > data2$co2_samp_2 ),18:20]) 
length(badBLrows)

data2[rownames(data2) %in% badBLrows,]

data2[which( data2$co2_samp_control > data2$co2_samp_1  | data2$co2_samp_control > data2$co2_samp_2 ),18] <- NA
 
data2$co2_samp_1_correct <- ifelse(is.na(data2$co2_samp_control), data2$co2_samp_1,data2$co2_samp_1 - data2$co2_samp_control)

data2$co2_samp_2_correct <- ifelse(is.na(data2$co2_samp_control), data2$co2_samp_2,data2$co2_samp_2 - data2$co2_samp_control)

data2$co2_final_samp <- pmax(data2$co2_samp_1_correct, data2$co2_samp_2_correct)

write.csv(data2, row.names = F, "data/data_final/corrected_ldeli_integrate.csv")

finaldata <- read.csv("data/data_final/corrected_ldeli_integrate.csv")

str(finaldata)

finaldata$date <- as.Date(finaldata$date, "%d/%m/%Y")
finaldata$lizmass <- finaldata$ch_lizmass - finaldata$ch_mass
finaldata$total_air <- finaldata$ch_vol - finaldata$lizmass
finaldata$total_co2 <- finaldata$frac_co2 * finaldata$total_air
finaldata$co2_pmin <- finaldata$total_co2 / finaldata$t_diff

write.csv(finaldata, row.names = F, "data/data_final/smr_rawfinal.csv")


