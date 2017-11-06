library(dplyr)
library(tidyr)

#Script to get indivual intercept and slope + their global intercept and slope and then plotting this as y
abline lines
#What do I need for this....
#Dataframe for each individual at each sampling period at each temp
#lizard.id, sampling.period, temp, pred.mr, intercept, slope
#so the intercept will be global.intercept at time + id.intercept blup at the time
#and the slope will be the global.slope + id.slope at the time - but I don't need a temp!? 


ab.predictions <- function(id = "ld0133", post, sampling.period=1){
  
  # Check that this combo is in the dataset
  focal.series <- paste(id, sampling.period, sep = "_")
  if(!(focal.series %in% dat$series)) return(NULL)
  print(focal.series)
  
  #setting up lizard and sampling period variables use for calculating predictions
  lizard.id <- paste("(Intercept).id.", id, sep = "")
  
  sampling.id <- paste("(Intercept).series.", id, "_", sampling.period, sep = "") # one can pick the sampling period, e.g. 1, or 2...10
  
  lizard.slope <- paste("inverseK_incb_temp.id.", id, sep = "")
  sampling.slope <- paste("inverseK_incb_temp.series.", id, "_", sampling.period, sep = "")
  
  post <- post %>% as.data.frame()
  
  #calculating predictions for my measured temperatures
  out <- do.call("cbind", lapply(1:length(incubation.temperatures), function(i){
    focal.answer <- post[, names(post) == "(Intercept)"] +            # global intercept
      post[, names(post) == lizard.id] + post[, names(post) == sampling.id] + # lizard- and sampling period-specific random intercepts
      (incubation.temperatures[i] * (post$inverseK_incb_temp + post[, names(post) == lizard.slope] + post[, names(post) == sampling.slope])) + # temp * (global.slope +  slope for THIS lizard in THIS sampling period)
      (post$inverseK_prior_temp2 * mean.prior.temp) # fixed effect of prior temp, assuming mean value of prior temp
    # we leave out mass, since average is zero, so we are calcualting this for an average mass liard
  }))
  
  #compiling predictions into data frame
  data.frame(pred.mr = c(out),
             Temperature = rep(incubation.temperatures, each = nrow(out)),
             Lizard = id,
             liz.intercept = post[, names(post) == "(Intercept)"] + post[, names(post) == lizard.id] + post[, names(post) == sampling.id] + (post$inverseK_prior_temp2 * mean.prior.temp),
             liz.slope = post$inverseK_incb_temp + post[, names(post) == lizard.slope] + post[, names(post) == sampling.slope],
             stringsAsFactors = FALSE)
}

# We will calculate the reaction norms assuming that all lizards had the mean value for 'prior temp'
mean.prior.temp <- mean(dat$inverseK_prior_temp2, na.rm = T)

# Here are the six incubator temperatures, 
incubation.temperatures <- unique(dat$inverseK_incb_temp)

#and the 42 lizard names
lizard.names <- data$id %>% as.character %>% unique %>% sort

m2 <- readRDS("output/rds/m2")

m2.S <- lapply(m2, function(m) m$Sol)
m2.Sol <- do.call(rbind, m2.S)
m2.S1 <- do.call(mcmc.list, m2.S)

m2.V <- lapply(m2, function(m) m$VCV)
m2.V1 <- do.call(mcmc.list, m2.V)
m2.VCV <- do.call(rbind, m2.V)

ab.output <- do.call("rbind", lapply(1:10, function(i){
  do.call("rbind", lapply(lizard.names, ab.predictions, post=cbind(m2.Sol, m2.VCV), sampling.period = i)) %>% mutate(sampling.period = i) %>% arrange(Temperature, pred.mr) %>% mutate(Lizard = factor(Lizard, levels = unique(Lizard)))
  }))

ab.reaction.norms <- ab.output %>% group_by(Temperature, Lizard, sampling.period) %>% summarise(pred.mr = posterior.mode(as.mcmc(pred.mr)), liz.intercept = posterior.mode(as.mcmc(liz.intercept)), liz.slope = posterior.mode(as.mcmc(liz.slope))) %>% as.data.frame
ab.reaction.norms$Temperature <- inverseK_to_C(ab.reaction.norms$Temperature)

str(ab.reaction.norms) ; head(ab.reaction.norms)

#ab.rn <- select(ab.reaction.norms, Temperature, Lizard, sampling.period, pred.mr) %>% mutate(series = paste0(Lizard, "_", sampling.period)) %>% arrange(series)
saveRDS(ab.rn, "output/rds/ab.rn")
ab.rn <- readRDS("output/rds/ab.rn")

lm_extract <- function(x){
    tmpmod <- lm(pred.mr ~ Temperature, data = x)
    coefs <-coef(tmpmod)
  return(coefs)
}

# split
  splitABRN <- split(ab.rn, ab.rn$series)

# apply
  ab.rn.coefs <- do.call(rbind, lapply(splitABRN, function(x) lm_extract(x)))
  ab.rn.coefs <- as.tbl(ab.rn.coefs)
  str(ab.rn.coefs)
  ab.rn.coefs$series <- rownames(ab.rn.coefs)
  ab.rn <- as.tbl(ab.rn)
  str(ab.rn)
  ab.rn <-left_join(ab.rn, ab.rn.coefs, by = "series")
  names(ab.rn)[1] <- "Temperature"
  names(ab.rn)[6] <- "b0"
  names(ab.rn)[7] <- "b1"

#pdf("output/fig/ab_reaction.norms.pdf", 10, 6)
ab.rn %>% 
  ggplot(aes(x = Temperature, y = pred.mr, group = Lizard)) + 
  geom_point(shape = 1, fill = "white", size = 1, color = "black") + 
  geom_abline(aes(intercept = b0, slope = b1, color = Lizard), alpha = 0.6) + 
  scale_x_continuous(breaks = c(22, 24, 26, 28, 30, 32)) + 
  facet_wrap(~ sampling.period, nrow = 2) + 
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  labs(x = expression(paste("Temperature ",degree,"C")), y = expression(Metabolic~rate~(CO[2]~min^{-1})))  
#dev.off()

#Careau code
plot(Phenotype.1~Phenotype.2,data=DATA,col=DATA$Individual,pch=1, xaxt="n",xlab="",ylab="",ylim=yrange,xlim=xrange)
for(i in unique(DATA$ID)){
  hold<-subset(DATA, ID==i)
  model<-lm(Phenotype.1~Phenotype.2,data=hold)
  clip(min(hold$Phenotype.2),max(hold$Phenotype.2),min(hold$Phenotype.1),max(hold$Phenotype.1))
  abline(model, col=i)
}
clip(min(DATA$Phenotype.2),max(DATA$Phenotype.2),min(DATA$Phenotype.1),max(DATA$Phenotype.1))
points(DATA.AVG$avg1~DATA.AVG$avg2, pch=16)
clip(min(DATA.AVG$avg2),max(DATA.AVG$avg2),min(DATA.AVG$avg1),max(DATA.AVG$avg1))
abline(lm(DATA.AVG$avg1~DATA.AVG$avg2), col="black",lwd=2)
mtext(expression(bold("A")),side=3, line=-1.5,adj=0.025)
mtext(expression(paste(italic("N")[ID], " = 100")), side=1,line=-1,adj=1)
mtext(expression(paste(italic("n")[trial], " = 20")), side=1,line=-2,adj=1)


#Graph to show between and within correltions between temps
lm_extract_mr <- function(x){
  tmpmod <- lm(mr.32 ~ mr.22, data = x)
  coefs <-coef(tmpmod)
  return(coefs)
}

#Filtering dataset to se Temps 32 and 22 only 
longfor.22.32 <- filter(ab.rn, Temperature == "32" | Temperature == "22") %>% spread(Temperature, pred.mr)
names(longfor.22.32)[6] <- "mr.22"
names(longfor.22.32)[7] <- "mr.32"

long.dat.22.32 <- select(longfor.22.32, Lizard, sampling.period, series, mr.22, mr.32)
#spliting dataset to calculating individual reg lines
split.22.32 <- split(long.dat.22.32, long.dat.22.32$Lizard)

dat.22.32.coefs <- as.data.frame(do.call("rbind", lapply(split.22.32, function(x) lm_extract_mr(x))))
dat.22.32.coefs$Lizard <- rownames(dat.22.32.coefs)
names(dat.22.32.coefs) <- c("b0", "b1", "Lizard")

str(dat.22.32.coefs) ; str(long.dat.22.32) #fucking dataframes, can't use left_join, convert as tbl
dat.22.32.coefs <- as.tbl(dat.22.32.coefs)
long.dat.22.32 <- as.tbl(long.dat.22.32)

long.dat.22.32.dat <- left_join(long.dat.22.32, dat.22.32.coefs, by = "Lizard")

ggplot(data = long.dat.22.32.dat, aes(y = mr.32, x = mr.22, group = "Lizard")) +
  geom_point(shape = 1, fill = "white", size = 1, color = "black") +
  #geom_abline(aes(intercept = b0, slope = b1, colour = Lizard)) +
  #stat_smooth(aes(group = Lizard, colour = Lizard), method = "lm", se = FALSE) + 
  geom_line(aes(group = Lizard, colour = Lizard), stat="smooth", method = "lm", alpha = 0.6) + 
  labs(y = expression(Metabolic~rate~(CO[2]~min^{-1})~at~32~paste(degree,"C")),
       x = expression(Metabolic~rate~(CO[2]~min^{-1})~at~22~paste(degree,"C"))) +
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())




ab.reaction.norms[ab.reaction.norms$Temperature == "32","pred.mr"]

