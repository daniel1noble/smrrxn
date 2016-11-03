setwd("~/Dropbox/smrrxn/")

maledat <- read.csv("data/LDMales_alive.csv")

#Constructing randomisation dataframe for male lizards

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


#Constructing randomisation dataframe for temperatures to incbuators

tempdat <- data.frame(period = rep(c(1:10), each = 24),
                      batch = rep(c(1:2), each = 12, times = 10),
                      day = rep(c(1,1,1,1,2,2,2,2,3,3,3,3), times = 10),
                      incubator = rep(c(1,2,1,2), times = 60))

tempdat <- tempdat[order(tempdat$incubator),]

write.csv(tempdat, file = "data/temprandom.csv", row.names = F)

#day 3 temp
replicate(20, sample(c(30,32)))

#day 2 temp
day_2 <- replicate(20, sample(c(26,28)))

#day 1 temp
day_1 <-replicate(20, sample(c(22,24)))

