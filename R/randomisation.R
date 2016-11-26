setwd("~/Dropbox/smrrxn/")

maledat <- read.csv("data/LDMales_alive.csv")
names(maledat)


#Constructing randomisation dataframe for male lizards

randat <- data.frame(period = rep(c(1:10), each = 42),
                     batch = rep(c(1:2), times = c(21, 21)),
                     incubator = rep(c(1,2,1,2), times = c(10,11,10,11)),
                     chamber = rep(c(1:21, 1:21), times = c(10)))

randat <- randat[order(randat$batch),]    
str(randat)

#Assigning males to B1 and B2

b1_ids <- maledat$id[1:21]
b2_ids <- maledat$id[22:42]

#Sampling B1 ids for 10 periods
allb1 <- replicate(10,sample(b1_ids, replace = F))
allb1 <- unlist(list(allb1))

allb2 <- replicate(10,sample(b2_ids, replace = F))
allb2 <- unlist(list(allb2))

randat$id <- c(allb1, allb2)

unique(randat[randat$batch == 1, 5]) 
unique(randat[randat$batch == 2, 5])

#merge randat and maledat together

randat$id <- as.factor(randat$id)

str(maledat)

final_randat <- merge(randat, maledat[,5:7])

final_randat<- final_randat[order(final_randat$period, final_randat$batch, final_randat$incubator, final_randat$chamber),]

final_randat <- final_randat[c(2,3,4,5,1,7,6)]

write.csv(final_randat, file = "data/randomisation_v3.csv", row.names = F) 


#Constructing randomisation dataframe for temperatures to incbuators

tempdat <- data.frame(period = rep(c(1:10), each = 24),
                      batch = rep(c(1:2), each = 12, times = 10),
                      day = rep(c(1,1,1,1,2,2,2,2,3,3,3,3), times = 10),
                      incubator = rep(c(1,2,1,2), times = 60))

tempdat <- tempdat[order(tempdat$incubator),]



#day 3 temp
replicate(20, sample(c(30,32)))

#day 2 temp
day_2 <- replicate(20, sample(c(26,28)))

#day 1 temp
day_1 <-replicate(20, sample(c(22,24)))

#randomise temps
tempdat$temp <- unlist(list(replicate(40, sample(c(22,24,26,28,30,32), replace =F))))

write.csv(tempdat, file = "data/temprandom_V2.csv", row.names = F)

getwd()
