setwd("~/Google Drive/PhD 2016 - 2020/Experiments/incubation_geneexp/")

inc1.eggsdat <- read.csv("incubator1_eggs_11.12.16.csv", header = T)
inc2.eggsdat <- read.csv("incubator2_eggs_11.12.16.csv", header = T)

str(inc1.eggsdat) ; names(inc1.eggsdat)
str(inc2.eggsdat) ; names(inc2.eggsdat)

eggscombined <- rbind(inc1.eggsdat, inc2.eggsdat)
str(eggscombined)

summary(inc1.eggsdat$egg_incub_days)
summary(inc2.eggsdat$egg_incub_days)

#Subsetting eggs from I1 that are less than 25 days old
inc1.eggsdat <- inc1.eggsdat[inc1.eggsdat$egg_incub_days < 25,]

#Subsetting eggs from I2 that are less than 39 days old
inc2.eggsdat <- inc2.eggsdat[inc2.eggsdat$egg_incub_days < 39,]

#sampling eggs for gene expression

#Inc 1
inc1.cull.eggs <- sample(inc1.eggsdat$egg_id,10, replace = F)
#egg0151 egg0149 egg0122 egg0141 egg0145 egg0093 egg0134 egg0152 egg0139 egg0130

inc2.cull.eggs <- sample(inc2.eggsdat$egg_id,10, replace = F)
#egg0044 egg0095 egg0084 egg0100 egg0121 egg0049 egg0089 egg0153 egg0132 egg0127

eggs_sample <- unlist(list(inc1.cull.eggs, inc2.cull.eggs))

eggculldat <- data.frame(egg_id = eggs_sample)

eggculldat_final <- merge(eggculldat, eggscombined)

#Checking

eggscombined[eggscombined$egg_id == "egg0132",]
eggscombined[eggscombined$egg_id == "egg0040",]

#What day to sample? for Inc 1

replicate(5, sample(seq(11, 19), 1, replace = F))
replicate(5, sample(seq(13, 22), 1, replace = F))

#What day to sample? for Inc 2

replicate(5, sample(seq(18, 32), 1, replace = F))
replicate(5, sample(seq(23, 36), 1, replace = F))

#Writing file

write.csv(eggculldat_final, file = "eggs_to_cull.csv", row.names = F)


