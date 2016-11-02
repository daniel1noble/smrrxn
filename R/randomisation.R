setwd("~/Dropbox/smrrxn/")

maledat <- read.csv("data/LDMales_alive.csv")

#Constructing randomisation dataframe

randat <- data.frame(period = rep(c(1:10), each = 45),
                     batch = rep(c(1:2), times = c(23, 22)),
                     incubator = rep(c(1,2,1,2), times = c(11,12,11,11)),
                     chamber = rep(c(1:23, 1:22), times = 10))
                              

#Assigning males to B1 and B2
ids <- replicate(10,sample(maledat$id, replace = F))
ids_all <- unlist(list(ids))

randat$id <- ids_all

table(randat[randat$period == "10",5])

#merge randat and maledat together
str(randat)
randat$id <- as.factor(randat$id)

str(maledat)

final_randat <- merge(randat, maledat[,1:3])

final_randat <- final_randat[c(2, 3,4,5,1,6,7)]

final_randat<- final_randat[order(final_randat$period, final_randat$batch, final_randat$incubator, final_randat$chamber),]

write.csv(final_randat, file = "data/randomisation.csv", row.names = F) 







######

B1_ids <- ids[1:23]
B2_ids <- ids[24:length(ids)]

#Assigning males to incubators for Batch 1

B1_incub_ids <- sample(B1_ids, replace = F)
B1_I1_ids <- B1_incub_ids[1:12]
B1_I2_ids <- B1_incub_ids[13:23]

B1_I1_chamber_ids <- sample(B1_I1_ids, replace = F)
B1_I2_chamber_ids <- sample(B1_I2_ids, replace = F)

B1_I1_chamber_ids <- as.character(B1_I1_chamber_ids)
B1_I2_chamber_ids  <- as.character(B1_I2_chamber_ids )

randat$male_id <- c(B1_I1_chamber_ids, B1_I2_chamber_ids)


#Assigning males to incubators in Batch 2

B2_incub_ids <- sample(B2_ids, replace = F)
B2_I1_ids <- incub_ids[1:12]
B2_I2_ids <- incub_ids[13:23]


