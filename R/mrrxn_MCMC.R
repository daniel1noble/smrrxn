setwd("~/Dropbox/smrrxn/")

#clear envir
rm(list = ls())

#load library
library(MCMCglmm)
library(parallel)
library(dplyr)
library(ggplot2)

source("R/functions/smr_functions.R")

#data

data <- read.csv("data/data_final/mrrxn_final_v2.csv")
data$id <- as.factor(data$id)
data$series <- as.factor(data$series)

varibs.need <- c("obs", "samp_period", "id" ,"batch", "series", "incb_num", "incb_temp_id", "defecate", "incb_temp", "z.incb_temp", "z.log.temp", "incb_temp_K", "inverseK_incb_temp", "body_temp", "z.body_temp", "z.log.body_temp", "body_temp_K", "inverseK_body_temp" , "z.prior_temp1", "z.log.prior_temp1", "prior_temp1_K", "inverseK_prior_temp1", "prior_temp2_K", "inverseK_prior_temp2", "z.prior_temp2", "z.log.prior_temp2", "orig_lizmass", "lizmass_nocombout", "log.mass", "z.log.mass","orig_co2_pmin", "co2pm_nocombout", "log.co2pmin", "z.log.co2pmin")

incl.vars <- names(data) %in% varibs.need
data <- data[incl.vars]

predictors <- c("z.log.co2pmin", "inverseK_incb_temp", "z.log.mass", "inverseK_prior_temp2", "id", "series")

dat <- data[complete.cases(data[,predictors]),]
str(dat)

# Define whether to run code

m1 <- FALSE
m2 <- FALSE
m3 <- FALSE
m4 <- FALSE

#priors
expanded.prior <- list(R = list(V = 1, nu = 0.002),
                       G = list(G1 = list(V = diag(2), nu = 0.002, alpha.V = diag(1000,2,2), alpha.mu = rep(0,2)),
                                G2 = list(V = diag(2), nu = 0.002, alpha.V = diag(1000,2,2), alpha.mu = rep(0,2))))

#m1
if(m1){
  m1 <- mclapply(1:3, function(i) {
    MCMCglmm(z.log.co2pmin ~ inverseK_incb_temp + z.log.mass + inverseK_prior_temp2,
             random = ~us(1+inverseK_incb_temp):id + us(1+inverseK_incb_temp):series,
             family = "gaussian",
             prior = expanded.prior,
             nitt = 7510000,
             burnin = 10000,
             thin = 5000,
             data = dat, 
             verbose = T)
  }, mc.cores = 3)
}else{
m1 <- readRDS("output/rds/m1")
}

#Calculating repeatabilities using m1
m1.S <- lapply(m1, function(m) m$Sol)
m1.Sol <- do.call(rbind, m1.S)
m1.S1 <- do.call(mcmc.list, m1.S)

m1.V <- lapply(m1, function(m) m$VCV)
m1.VCV <- do.call(rbind, m1.V)
m1.V1 <- do.call(mcmc.list, m1.V)

gelman.diag(m1.S1)
summary(m1.Sol)
posterior.mode(m1.Sol)
HPDinterval(as.mcmc(m1.Sol))

gelman.diag(m1.V1, multivariate = F)
summary(m1.VCV)
posterior.mode(m1.VCV)
HPDinterval(as.mcmc(m1.VCV))

#Tabulating the model output using m1
Table1 <- data.frame(matrix(nrow = 15 , ncol = 3))
#rownames(Table1) <- c("Interpcet", "inverseK_incb_temp", "z.log.mass", "inverseK_prior_temp2",
#                      "IDintercept", "IDslope", "COV IDintercept-IDslope",
#                      "Sintercept", "Sslope", "COV Sintercept-Sslope",
#                      "e", "Rint", "Rslope", "Rshort", "Rlong")

rownames(Table1) <- c(names(c(posterior.mode(m1.Sol), posterior.mode(m1.VCV)[c(1,4,2,5,8,6,9)])), "Rint", "Rslope", "Rshort", "Rlong")
colnames(Table1) <- c("estimate", "lower", "upper")

#Tabulating fixed efs and CIs
Table1[1:4, 1] <- posterior.mode(m1.Sol)
Table1[1:4, 2:3] <- HPDinterval(as.mcmc(m1.Sol))

#Tabulating ran efs and CIs
Table1[5:11, 1] <- posterior.mode(m1.VCV)[c(1,4,2,5,8,6,9)]
Table1[5:11, 2:3] <- HPDinterval(as.mcmc(m1.VCV))[c(1,4,2,5,8,6,9),]

#Tabulating repeatabilities and CIs
head(m1.VCV)
colnames(m1.VCV)

#Rint
R.int <- m1.VCV[,"(Intercept):(Intercept).id"] / ( m1.VCV[,"(Intercept):(Intercept).id"] + m1.VCV[,"(Intercept):(Intercept).series"] ) 

Table1[12,1] <- posterior.mode(R.int) #R.int
Table1[12,2:3] <- HPDinterval(as.mcmc(R.int)) #R.int CIs

#Rslope
R.slope <- m1.VCV[,"inverseK_incb_temp:inverseK_incb_temp.id"] / ( m1.VCV[,"inverseK_incb_temp:inverseK_incb_temp.id"] + m1.VCV[,"inverseK_incb_temp:inverseK_incb_temp.series"] )

Table1[13,1] <-  posterior.mode(R.slope) #Rslope
Table1[13,2:3] <- HPDinterval(as.mcmc(R.slope)) #Rslope CIs

#Rshort in intercept
R.short <- ( m1.VCV[,"(Intercept):(Intercept).id"] + m1.VCV[,"(Intercept):(Intercept).series"] ) / ( m1.VCV[,"(Intercept):(Intercept).id"] + m1.VCV[,"(Intercept):(Intercept).series"]  + m1.VCV[,"units"] )

Table1[14,1] <- posterior.mode(R.short) #R.short
Table1[14,2:3] <- HPDinterval(as.mcmc(R.short)) #R.short CIs

#Rlong in intercept
R.long <-  m1.VCV[,"(Intercept):(Intercept).id"] / ( m1.VCV[,"(Intercept):(Intercept).id"] + m1.VCV[,"(Intercept):(Intercept).series"]  + m1.VCV[,"units"] )

Table1[15,1] <-  posterior.mode(R.long) #Rlong
Table1[15,2:3] <- HPDinterval(as.mcmc(R.long)) #Rlong CIs

write.csv(Table1, "output/table/Table1.csv")
write.csv(round(Table1,2), "output/table/Table1_rounded.csv")

#Calculating correlation of ID:slope and series:slope
id.slope.cor <- m1.VCV[,"inverseK_incb_temp:(Intercept).id"] / ( sqrt(m1.VCV[,"(Intercept):(Intercept).id"]) * sqrt(m1.VCV[,"inverseK_incb_temp:inverseK_incb_temp.id"]) )
hist(id.slope.cor)
posterior.mode(id.slope.cor)
HPDinterval(as.mcmc(id.slope.cor))

series.slope.cor <- m1.VCV[,"inverseK_incb_temp:(Intercept).series"] / ( sqrt(m1.VCV[,"(Intercept):(Intercept).series"]) * sqrt(m1.VCV[,"inverseK_incb_temp:inverseK_incb_temp.series"]) )
hist(series.slope.cor)
posterior.mode(series.slope.cor)
HPDinterval(as.mcmc(series.slope.cor))

#m2
if(m2){
  m2 <- mclapply(1:3, function(i) {
    MCMCglmm(z.log.co2pmin ~ inverseK_incb_temp + z.log.mass + inverseK_prior_temp2,
             random = ~us(1+inverseK_incb_temp):id + us(1+inverseK_incb_temp):series,
             family = "gaussian",
             prior = expanded.prior,
             nitt = 7510000,
             burnin = 10000,
             thin = 5000,
             pr = TRUE,
             data = dat, 
             verbose = T)
  }, mc.cores = 3)
}else{
  m2 <- readRDS("output/rds/m2")
}

m2.S <- lapply(m2, function(m) m$Sol)
m2.Sol <- do.call(rbind, m2.S)
m2.S1 <- do.call(mcmc.list, m2.S)

gelman.diag(m2.S1, multivariate = F)
summary(m2.Sol)
posterior.mode(m2.Sol)
HPDinterval(as.mcmc(m2.Sol))

m2.V <- lapply(m2, function(m) m$VCV)
m2.V1 <- do.call(mcmc.list, m2.V)
m2.VCV <- do.call(rbind, m2.V)

plot(m2.VCV)
gelman.diag(m2.VCV, multivariate = F)
summary(m2.VCV)
posterior.mode(m2.VCV)
HPDinterval(as.mcmc(m2.VCV))

#Individual predictions from m2 for plots
get.predictions <- function(id = "ld0133", post, sampling.period=1){
  
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
  data.frame(predicted = c(out),
             Temperature = rep(incubation.temperatures, each = nrow(out)),
             Lizard = id, 
             stringsAsFactors = FALSE)
}

# We will calculate the reaction norms assuming that all lizards had the mean value for 'prior temp'
mean.prior.temp <- mean(dat$inverseK_prior_temp2, na.rm = T)

# Here are the six incubator temperatures, 
incubation.temperatures <- unique(dat$inverseK_incb_temp)

#and the 42 lizard names
lizard.names <- data$id %>% as.character %>% unique %>% sort

output <- do.call("rbind", lapply(1:10, function(i){
  do.call("rbind", lapply(lizard.names, get.predictions, post=cbind(m2.Sol, m2.VCV), sampling.period = i)) %>% mutate(sampling.period = i) %>% arrange(Temperature, predicted) %>% mutate(Lizard = factor(Lizard, levels = unique(Lizard))) 
}))
#saveRDS(output, "output/id.rxnnorm.preds")

output <- readRDS("output/id.rxnnorm.preds")
reaction.norms <- output %>% group_by(Temperature, Lizard, sampling.period) %>% summarise(posterior.mode = posterior.mode(as.mcmc(predicted)), lowerCI = as.numeric(HPDinterval(as.mcmc(predicted)))[1], upperCI = as.numeric(HPDinterval(as.mcmc(predicted)))[2])

reaction.norms$Temperature <- inverseK_to_C(reaction.norms$Temperature)
str(reaction.norms)

#Plotting reactions
#reaction.norms %>% ggplot(aes(Lizard, posterior.mode)) + geom_hline(yintercept = 0, linetype = 2) + geom_errorbar(aes(ymin = lowerCI, ymax = upperCI), width = 0) + geom_point() + coord_flip() + facet_wrap(~Temperature) + ylab("Respiration rate")

#pdf("output/fig/reaction.norms.pdf", 10, 6)
reaction.norms %>% 
  ggplot(aes(x = Temperature, y = posterior.mode, group = Lizard)) +
  geom_line(color = "grey") + 
  geom_point(shape = 1, fill = "white", size = 1, color = "black") +
  facet_wrap(~ sampling.period, nrow = 2) + 
  scale_x_continuous(breaks = c(22, 24, 26, 28, 30, 32)) + 
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  labs(x = expression(paste("Temperature ",degree,"C")), y = expression(Metabolic~rate~(CO[2]~min^{-1})))
#dev.off()

#pdf("output/fig/ID.rxn.norm.pdf", 9,9)
reaction.norms %>% 
  ggplot(aes(x = Temperature, y = posterior.mode, group = sampling.period, color = sampling.period)) + 
  geom_line() + 
  geom_point(shape = 1, fill = "white", size = 1, color = "black") + 
  scale_x_continuous(breaks = c(22, 24, 26, 28, 30, 32)) + 
  facet_wrap(~ Lizard) + 
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  labs(x = expression(paste("Temperature ",degree,"C")), y = expression(Metabolic~rate~(CO[2]~min^{-1})))
#dev.off()

#Plotting covariances
#Individual predictions from m2 for plots
get.intercept.slope <- function(id = "ld0133", post, sampling.period=1){
  
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
  
  #compiling into data frame
  out <- data.frame(Lizard = id,
                    global.int = post[, names(post) == "(Intercept)"], # global intercept
                    ind.int = post[, names(post) == lizard.id], #lizard- random intercepts
                    series.int = post[, names(post) == sampling.id], # sampling period-specific random intercepts
                    ind.slope = post[, names(post) == lizard.slope],
                    series.slope = post[, names(post) == sampling.slope],
                    stringsAsFactors = FALSE)
}


#and the 42 lizard names
lizard.names <- dat$id %>% as.character %>% unique %>% sort

#covar.output <- do.call("rbind", lapply(1:10, function(i){
  #do.call("rbind", lapply(lizard.names, get.intercept.slope, post=cbind(m2.Sol, m2.VCV), sampling.period = i)) %>% mutate(sampling.period = i) %>% arrange(sampling.period) %>% mutate(Lizard = factor(Lizard, levels = unique(Lizard))) 
#}))
#saveRDS(covar.output, "output/covar.preds")

covar.output <- readRDS("output/covar.preds")
cor.int.slop <- covar.output %>% group_by(Lizard, sampling.period) %>% summarise(global.int = posterior.mode(as.mcmc(global.int)),
                                                                           Ind.int = posterior.mode(as.mcmc(ind.int)),
                                                                           Series.int = posterior.mode(as.mcmc(series.int)),
                                                                           Ind.slope = posterior.mode(as.mcmc(ind.slope)),
                                                                           Series.slope = posterior.mode(as.mcmc(series.slope)),
                                                                           lower.Ind.int = as.numeric(HPDinterval(as.mcmc(ind.int))[,1]),
                                                                           upper.Ind.int = as.numeric(HPDinterval(as.mcmc(ind.int))[,2]),
                                                                           lower.Series.int = HPDinterval(as.mcmc(series.int))[,1],
                                                                           upper.Series.int = HPDinterval(as.mcmc(series.int))[,2],
                                                                           lower.Ind.slope = HPDinterval(as.mcmc(ind.slope))[,1],
                                                                           upper.Ind.slope = HPDinterval(as.mcmc(ind.slope))[,2],
                                                                           lower.Series.slope = HPDinterval(as.mcmc(series.slope))[,1],
                                                                           upper.Series.slope = HPDinterval(as.mcmc(series.slope))[,2],)

series.mass <- select(dat, id, samp_period, orig_lizmass) %>% 
  group_by(id, samp_period) %>% 
  summarise(mean_mass = mean(orig_lizmass))

series.mass$z.log.mass <- scale(log(series.mass$mean_mass))
names(series.mass) <- c("Lizard", "sampling.period", "mean_mass", "z.log.mass" )
names(cor.int.slop)

cor.int.slop2 <- left_join(cor.int.slop, series.mass)
str(cor.int.slop2)



#Plotting ints with slopes

#pdf("output/fig/covariance.ID.int.slope.pdf", 10, 6)
#fig2a <- 
  cor.int.slop %>% ggplot(aes(y = Ind.int, x = Ind.slope)) + 
  geom_errorbar(aes(x = Ind.slope, ymin = lower.Ind.int, ymax = upper.Ind.int), color = "#999999", width = 0, size = 0.25) + 
  geom_errorbarh(aes(x = Ind.slope, xmin = lower.Ind.slope, xmax = upper.Ind.slope), color = "#999999", size = 0.25) + 
  geom_point(shape = 1, fill = "white", size = 1, color = "black")  + facet_wrap(~ sampling.period, nrow = 2) + 
  theme_bw() + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = "Slope for ID", y = "Intercept for ID") 
#dev.off()
  
  #pdf("output/fig/covariance.ID.int.slope.pdf", 10, 6)
  #fig2a <- 
  cor.int.slop2 %>% ggplot(aes(y = Ind.int, x = Ind.slope, color = z.log.mass)) + 
    geom_errorbar(aes(x = Ind.slope, ymin = lower.Ind.int, ymax = upper.Ind.int), width = 0, size = 0.25) + 
    geom_errorbarh(aes(x = Ind.slope, xmin = lower.Ind.slope, xmax = upper.Ind.slope), size = 0.25) + 
    geom_point(shape = 1, size = 1)  + 
    facet_wrap(~ sampling.period, nrow = 2) + 
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    labs(x = "Slope for ID", y = "Intercept for ID") 
  #dev.off()

#pdf("output/fig/covariance.series.int.slope.pdf", 10, 6)
#fig2b <- 
  cor.int.slop %>% ggplot(aes(y = Series.int, x = Series.slope)) + 
  geom_errorbar(aes(x = Series.slope, ymin = lower.Series.int, ymax = upper.Series.int), width = 0, size = 0.25, color = "#999999") +
  geom_errorbarh(aes(x = Series.slope, xmin = lower.Series.slope, xmax = upper.Series.slope),  color = "#999999", size = 0.25) + 
  geom_point(shape = 1, fill = "white", size = 1, color = "black")  + 
  facet_wrap(~ sampling.period, nrow = 2) +  
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = "Slope for Series", y = "Intercept for Series") 
#dev.off()

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#pdf("output/fig/covariance.ID.Series.pdf", 10, 10)
multiplot(fig2a, fig2b)
#dev.off()

#pdf("output/fig/covariance.series.int.slope.samp.period.pdf", 10, 10)
cor.int.slop2 %>% 
  ggplot(aes(y = Series.int, x = Series.slope, color = sampling.period)) + 
  geom_errorbar(aes(x = Series.slope, ymin = lower.Series.int, ymax = upper.Series.int), width = 0, size = 0.8) + 
  geom_errorbarh(aes(x = Series.slope, xmin = lower.Series.slope, xmax = upper.Series.slope), size = 0.8) + 
  geom_point(shape = 1, fill = "white", size = 1, color = "black")  + 
  facet_wrap(~ Lizard) +  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = "Slope for Series", y = "Intercept for Series") 
#dev.off()


#pdf("output/fig/", 10, 10)
cor.int.slop2 %>% 
  ggplot(aes(y = Series.int, x = Series.slope, color = z.log.mass)) + 
  geom_errorbar(aes(x = Series.slope, ymin = lower.Series.int, ymax = upper.Series.int), width = 0, size = 0.8) + 
  geom_errorbarh(aes(x = Series.slope, xmin = lower.Series.slope, xmax = upper.Series.slope), size = 0.8) + 
  geom_point(shape = 1, fill = "white", size = 1, color = "black")  + 
  facet_wrap(~ Lizard) +  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = "Slope for Series", y = "Intercept for Series") +
  scale_colour_gradient(low = "paleturquoise", high = "paleturquoise4", guide = "colourbar")
#dev.off()

#plot series changes in intercept with mass for each lizard

fig3a <- cor.int.slop2 %>% 
  ggplot(aes(y = Series.int, x = z.log.mass, color = sampling.period)) + 
  geom_errorbar(aes(x = z.log.mass, ymin = lower.Series.int, ymax = upper.Series.int), width = 0, size = 0.3) +
  #geom_point(shape = 1, size = 1, color = "black", fill = "white")  + 
  geom_point(shape = 1, size = 1)  + 
  facet_wrap(~ Lizard) +  
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = "", y = "Intercept for series") 
  #labs(x = "Z-transformed log mass (g)", y = "Intercept for series") 

fig3b <- cor.int.slop2 %>% 
  ggplot(aes(y = Series.slope, x = z.log.mass, color = sampling.period)) + 
  geom_errorbar(aes(x = z.log.mass, ymin = lower.Series.slope, ymax = upper.Series.slope), size = 0.3) + 
  #geom_point(shape = 1, size = 1, color = "black", fill = "white")  + 
  geom_point(shape = 1, size = 1)  + 
  facet_wrap(~ Lizard) +  
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = "Z-transformed log mass (g)", y = "Slope for series") 

#pdf("output/fig/covariance.ID.Series.Mass.pdf", 12, 20)
multiplot(fig3a, fig3b)
#dev.off()

#Plotting across sampling period

cor.int.slop2 %>% 
  ggplot(aes(y = Series.int, x = z.log.mass)) + 
  geom_errorbar(aes(x = z.log.mass, ymin = lower.Series.int, ymax = upper.Series.int), width = 0, size = 0.3) + 
  #geom_point(shape = 1, size = 1, color = "black", fill = "white")  + 
  geom_point(shape = 1, size = 1)  + 
  facet_wrap(~sampling.period, nrow = 2) +  
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = "Z-transformed log mass (g)", y = "Intercept for series") 

cor.int.slop2 %>% 
  ggplot(aes(y = Series.slope, x = z.log.mass)) + 
  geom_errorbar(aes(x = z.log.mass, ymin = lower.Series.slope, ymax = upper.Series.slope), size = 0.3) + 
  #geom_point(shape = 1, size = 1, color = "black", fill = "white")  + 
  geom_point(shape = 1, size = 1)  + 
  facet_wrap(~sampling.period, nrow = 2) +  
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = "Z-transformed log mass (g)", y = "Slope for series") 

#PLotting ID intercepts and slopes with mass 
#pdf("output/fig/ID.covar.mass.pdf", 10, 6)
cor.int.slop2 %>% ggplot(aes(y = Ind.int, x = Ind.slope, color = z.log.mass)) + 
  geom_errorbar(aes(x = Ind.slope, ymin = lower.Ind.int, ymax = upper.Ind.int), width = 0, size = 0.25) + 
  geom_errorbarh(aes(x = Ind.slope, xmin = lower.Ind.slope, xmax = upper.Ind.slope), size = 0.25) + 
  geom_point(shape = 1, size = 1)  + 
  facet_wrap(~ sampling.period, nrow = 2) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = "Slope for ID", y = "Intercept for ID")
#dev.off()

fig4a <- cor.int.slop2 %>% ggplot(aes(y = Ind.int, x = z.log.mass, color = Lizard)) + 
  geom_errorbar(aes(x = z.log.mass, ymin = lower.Ind.int, ymax = upper.Ind.int), width = 0, size = 0.25) + 
  geom_point(shape = 1, size = 1)  + 
  facet_wrap(~ sampling.period, nrow = 2) + 
  theme_bw() + 
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = "Z-transformed log mass (g)", y = "Intercept for ID")

fig4b <- cor.int.slop2 %>% ggplot(aes(y = Ind.slope, x = z.log.mass, color = Lizard)) + 
  geom_errorbar(aes(x = z.log.mass, ymin = lower.Ind.slope, ymax = upper.Ind.slope), size = 0.25) + 
  geom_point(shape = 1, size = 1)  + 
  facet_wrap(~ sampling.period, nrow = 2) + 
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = "Z-transformed log mass (g)", y = "Slope for ID") 

#pdf("output/fig/between.ID.mass.int.slop.pdf", 10, 10)
multiplot(fig4a, fig4b)
#dev.off()
#Temperature specific repeatability 
#m4

multi.prior <- list(R = list(V = diag(6), nu = 0.01), G = list(G1 = list(V = diag(6), nu = 0.01)))

if(m4){
m4 <- mclapply(1:3, function(i) {
  MCMCglmm(cbind(t_22, t_24, t_26, t_28, t_30, t_32) ~  z.log.mass,
           random= ~us(trait):id,
           rcov = ~us(trait):units,
           family = c(rep("gaussian", 6)),
           prior = multi.prior,
           nitt = 7510000,
           burnin = 10000,
           thin = 5000,
           data = data, 
           verbose = T)
}, mc.cores = 3)
} else{
  m4 <- readRDS("output/rds/m4_v2")
}

m4.S <- lapply(m4, function(m) m$Sol)
m4.Sol <- do.call(rbind, m4.S)
m4.S1 <- do.call(mcmc.list, m4.S)

gelman.diag(m4.S1, multivariate = F)
plot(m4.S1)
summary(m4.S1)
posterior.mode(m4.Sol)
HPDinterval(as.mcmc(m4.Sol))

m4.V <- lapply(m4, function(m) m$VCV)
m4.V1 <- do.call(mcmc.list, m4.V)
m4.VCV <- do.call(rbind, m4.V)

gelman.diag(m4.V1, multivariate = F)
summary(m4.V1)
posterior.mode(m4.VCV)
HPDinterval(as.mcmc(m4.VCV))

#Tabulating the model output using m1
Table2 <- data.frame(matrix(nrow = 6 , ncol = 6))
rownames(Table2) <- c(sort(unique(dat$incb_temp)))
colnames(Table2) <- c(sort(unique(dat$incb_temp)))

#Variances - the diagonal of matrix
names(posterior.mode(m4.VCV))

Table2[1,1] <- posterior.mode(m4.VCV)["traitt_22:traitt_22.id"]
Table2[2,2] <- posterior.mode(m4.VCV)["traitt_24:traitt_24.id"]
Table2[3,3] <- posterior.mode(m4.VCV)["traitt_26:traitt_26.id"]
Table2[4,4] <- posterior.mode(m4.VCV)["traitt_28:traitt_28.id"]
Table2[5,5] <- posterior.mode(m4.VCV)["traitt_30:traitt_30.id"]
Table2[6,6] <- posterior.mode(m4.VCV)["traitt_32:traitt_32.id"]

#Covariance - upper off diagonal of matrix 
#row of 22 
Table2[1,2] <- posterior.mode(m4.VCV)["traitt_22:traitt_24.id"]
Table2[1,3] <- posterior.mode(m4.VCV)["traitt_22:traitt_26.id"]
Table2[1,4] <- posterior.mode(m4.VCV)["traitt_22:traitt_28.id"]
Table2[1,5] <- posterior.mode(m4.VCV)["traitt_22:traitt_30.id"]
Table2[1,6] <- posterior.mode(m4.VCV)["traitt_22:traitt_32.id"]

#row of 24 
Table2[2,3] <- posterior.mode(m4.VCV)["traitt_24:traitt_26.id"]
Table2[2,4] <- posterior.mode(m4.VCV)["traitt_24:traitt_28.id"]
Table2[2,5] <- posterior.mode(m4.VCV)["traitt_24:traitt_30.id"]
Table2[2,6] <- posterior.mode(m4.VCV)["traitt_24:traitt_32.id"]

#row of 26 
Table2[3,4] <- posterior.mode(m4.VCV)["traitt_26:traitt_28.id"]
Table2[3,5] <- posterior.mode(m4.VCV)["traitt_26:traitt_30.id"]
Table2[3,6] <- posterior.mode(m4.VCV)["traitt_26:traitt_32.id"]

#row of 28 
Table2[4,5] <- posterior.mode(m4.VCV)["traitt_28:traitt_30.id"]
Table2[4,6] <- posterior.mode(m4.VCV)["traitt_28:traitt_32.id"]

#row of 30 
Table2[5,6] <- posterior.mode(m4.VCV)["traitt_30:traitt_32.id"]

#Correlation - lower off diagonal of matrix 
#col of 22
Table2[2,1] <- posterior.mode(m4.VCV[,"traitt_22:traitt_24.id"] / (sqrt(m4.VCV[,"traitt_22:traitt_22.id"]) * sqrt(m4.VCV[,"traitt_24:traitt_24.id"]))) 
Table2[3,1] <- posterior.mode(m4.VCV[,"traitt_22:traitt_26.id"] / (sqrt(m4.VCV[,"traitt_22:traitt_22.id"]) * sqrt(m4.VCV[,"traitt_26:traitt_26.id"]))) 
Table2[4,1] <- posterior.mode(m4.VCV[,"traitt_22:traitt_28.id"] / (sqrt(m4.VCV[,"traitt_22:traitt_22.id"]) * sqrt(m4.VCV[,"traitt_28:traitt_28.id"]))) 
Table2[5,1] <- posterior.mode(m4.VCV[,"traitt_22:traitt_30.id"] / (sqrt(m4.VCV[,"traitt_22:traitt_22.id"]) * sqrt(m4.VCV[,"traitt_30:traitt_30.id"])))
Table2[6,1] <- posterior.mode(m4.VCV[,"traitt_22:traitt_32.id"] / (sqrt(m4.VCV[,"traitt_22:traitt_22.id"]) * sqrt(m4.VCV[,"traitt_32:traitt_32.id"]))) 

#col of 24
Table2[3,2] <- posterior.mode(m4.VCV[,"traitt_24:traitt_26.id"] / (sqrt(m4.VCV[,"traitt_24:traitt_24.id"]) * sqrt(m4.VCV[,"traitt_26:traitt_26.id"]))) 
Table2[4,2] <- posterior.mode(m4.VCV[,"traitt_24:traitt_28.id"] / (sqrt(m4.VCV[,"traitt_24:traitt_24.id"]) * sqrt(m4.VCV[,"traitt_28:traitt_28.id"]))) 
Table2[5,2] <- posterior.mode(m4.VCV[,"traitt_24:traitt_30.id"] / (sqrt(m4.VCV[,"traitt_24:traitt_24.id"]) * sqrt(m4.VCV[,"traitt_30:traitt_30.id"])))
Table2[6,2] <- posterior.mode(m4.VCV[,"traitt_24:traitt_32.id"] / (sqrt(m4.VCV[,"traitt_24:traitt_24.id"]) * sqrt(m4.VCV[,"traitt_32:traitt_32.id"]))) 

#col of 26
Table2[4,3] <- posterior.mode(m4.VCV[,"traitt_26:traitt_28.id"] / (sqrt(m4.VCV[,"traitt_26:traitt_26.id"]) * sqrt(m4.VCV[,"traitt_28:traitt_28.id"]))) 
Table2[5,3] <- posterior.mode(m4.VCV[,"traitt_26:traitt_30.id"] / (sqrt(m4.VCV[,"traitt_26:traitt_26.id"]) * sqrt(m4.VCV[,"traitt_30:traitt_30.id"])))
Table2[6,3] <- posterior.mode(m4.VCV[,"traitt_26:traitt_32.id"] / (sqrt(m4.VCV[,"traitt_26:traitt_26.id"]) * sqrt(m4.VCV[,"traitt_32:traitt_32.id"]))) 

#col of 28
Table2[5,4] <- posterior.mode(m4.VCV[,"traitt_28:traitt_30.id"] / (sqrt(m4.VCV[,"traitt_28:traitt_28.id"]) * sqrt(m4.VCV[,"traitt_30:traitt_30.id"])))
Table2[6,4] <- posterior.mode(m4.VCV[,"traitt_28:traitt_32.id"] / (sqrt(m4.VCV[,"traitt_28:traitt_28.id"]) * sqrt(m4.VCV[,"traitt_32:traitt_32.id"])))

#col of 30
Table2[6,5] <- posterior.mode(m4.VCV[,"traitt_30:traitt_32.id"] / (sqrt(m4.VCV[,"traitt_30:traitt_30.id"]) * sqrt(m4.VCV[,"traitt_32:traitt_32.id"]))) 

write.csv(Table2, row.names = F, "output/table/Table2.csv")
write.csv(round(Table2, 2), row.names = F, "output/table/Table2_rounded.csv")

#Temperature specfic repeatability
Table3 <- data.frame(matrix(nrow = 6 , ncol = 3))
rownames(Table3) <- c(sort(unique(dat$incb_temp)))
colnames(Table3) <- c("estimate", "lower", "upper")

#22
Table3[1,1] <- posterior.mode(m4.VCV[,"traitt_22:traitt_22.id"] / ( m4.VCV[,"traitt_22:traitt_22.id"] + m4.VCV[,"traitt_22:traitt_22.units"]))
Table3[1,2:3] <- HPDinterval(as.mcmc(m4.VCV[,"traitt_22:traitt_22.id"] / ( m4.VCV[,"traitt_22:traitt_22.id"] + m4.VCV[,"traitt_22:traitt_22.units"])))

#24
Table3[2,1] <- posterior.mode(m4.VCV[,"traitt_24:traitt_24.id"] / ( m4.VCV[,"traitt_24:traitt_24.id"] + m4.VCV[,"traitt_24:traitt_24.units"]))
Table3[2,2:3] <- HPDinterval(as.mcmc(m4.VCV[,"traitt_24:traitt_24.id"] / ( m4.VCV[,"traitt_24:traitt_24.id"] + m4.VCV[,"traitt_24:traitt_24.units"])))

#26
Table3[3,1] <- posterior.mode(m4.VCV[,"traitt_26:traitt_26.id"] / ( m4.VCV[,"traitt_26:traitt_26.id"] + m4.VCV[,"traitt_26:traitt_26.units"]))
Table3[3,2:3] <- HPDinterval(as.mcmc(m4.VCV[,"traitt_26:traitt_26.id"] / ( m4.VCV[,"traitt_26:traitt_26.id"] + m4.VCV[,"traitt_26:traitt_26.units"])))

#28
Table3[4,1] <- posterior.mode(m4.VCV[,"traitt_28:traitt_28.id"] / ( m4.VCV[,"traitt_28:traitt_28.id"] + m4.VCV[,"traitt_28:traitt_28.units"]))
Table3[4,2:3] <- HPDinterval(as.mcmc(m4.VCV[,"traitt_28:traitt_28.id"] / ( m4.VCV[,"traitt_28:traitt_28.id"] + m4.VCV[,"traitt_28:traitt_28.units"])))

#30
Table3[5,1] <- posterior.mode(m4.VCV[,"traitt_30:traitt_30.id"] / ( m4.VCV[,"traitt_30:traitt_30.id"] + m4.VCV[,"traitt_30:traitt_30.units"]))
Table3[5,2:3] <- HPDinterval(as.mcmc(m4.VCV[,"traitt_30:traitt_30.id"] / ( m4.VCV[,"traitt_30:traitt_30.id"] + m4.VCV[,"traitt_30:traitt_30.units"])))

#32
Table3[6,1] <- posterior.mode(m4.VCV[,"traitt_32:traitt_32.id"] / ( m4.VCV[,"traitt_32:traitt_32.id"] + m4.VCV[,"traitt_32:traitt_32.units"]))
Table3[6,2:3] <- HPDinterval(as.mcmc(m4.VCV[,"traitt_32:traitt_32.id"] / ( m4.VCV[,"traitt_32:traitt_32.id"] + m4.VCV[,"traitt_32:traitt_32.units"])))

write.csv(Table3, row.names = F, 'output/table/Table3.csv')
write.csv(round(Table3, 2), row.names = F, 'output/table/Table3_rounded.csv')

#Plotting thermal repeatability
table3.dat <- Table3
table3.dat$temp <- row.names(Table3)
colnames(table3.dat)[1] <- "Repeatability"

#pdf("output/fig/thermal.rep.pdf", 6, 6)
table3.dat %>% ggplot(aes(y = Repeatability, x = temp)) + 
  geom_errorbar(aes(x = temp, ymin = lower, ymax = upper), size = 0.4, width = 0.2) + 
  geom_point(shape = 1, size = 1)  + 
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = expression(paste("Temperature ",degree,"C")), y = "Repeatability") 
#dev.off()
