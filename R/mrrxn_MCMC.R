setwd("~/Dropbox/smrrxn/")

#clear envir
rm(list = ls())

#load library
library(MCMCglmm)
library(parallel)
library(dplyr)
library(ggplot2)
library(reshape2)

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
m1 <- readRDS("output/rds/m1.log")
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

#output <- do.call("rbind", lapply(1:10, function(i){
  #do.call("rbind", lapply(lizard.names, get.predictions, post=cbind(m2.Sol, m2.VCV), sampling.period = i)) %>% mutate(sampling.period = i) %>% arrange(Temperature, predicted) %>% mutate(Lizard = factor(Lizard, levels = unique(Lizard)))
#}))

#saveRDS(output, "output/id.rxnnorm.preds")
output <- readRDS("output/id.rxnnorm.preds")

reaction.norms <- output %>% group_by(Temperature, Lizard, sampling.period) %>% summarise(posterior.mode = posterior.mode(as.mcmc(predicted)), lowerCI = as.numeric(HPDinterval(as.mcmc(predicted)))[1], upperCI = as.numeric(HPDinterval(as.mcmc(predicted)))[2])

reaction.norms$Temperature <- inverseK_to_C(reaction.norms$Temperature)
str(reaction.norms)

#Plotting reactions
#reaction.norms %>% ggplot(aes(Lizard, posterior.mode)) + geom_hline(yintercept = 0, linetype = 2) + geom_errorbar(aes(ymin = lowerCI, ymax = upperCI), width = 0) + geom_point() + coord_flip() + facet_wrap(~Temperature) + ylab("Respiration rate")

#pdf("output/fig/reaction.norms.pdf", 10, 6)
reaction.norms %>% 
  ggplot(aes(x = Temperature, y = posterior.mode, group = Lizard, color = Lizard)) +
  geom_line() + 
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
pdf("output/fig/AES_Fig1.pdf", 9,9)
filter(reaction.norms, sampling.period ==5) %>% 
  ggplot(aes(x = Temperature, y = posterior.mode, group = Lizard, color = Lizard)) +
  geom_line(aes(group = Lizard, colour = Lizard), stat="smooth", method = "lm", alpha = 0.6, lwd = 1.5) +
  geom_point(shape = 1, size = 1, colour = "black") +
  #facet_wrap(~ sampling.period, nrow = 2) + 
  scale_x_continuous(breaks = c(22, 24, 26, 28, 30, 32)) + 
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold")) + 
  labs(x = expression(paste("Temperature ",degree,"C")), y = expression(z-transformed~Metabolic~rate~(VCO[2]~min^{-1})))
dev.off()

pdf("output/fig/AES_Fig2.pdf", 9,9)
filter(reaction.norms, sampling.period ==5) %>% 
  ggplot(aes(x = Temperature, y = posterior.mode, group = Lizard, color = Lizard)) +
  #geom_line(aes(group = Lizard, colour = Lizard), stat="smooth", method = "lm", alpha = 0.6, lwd = 1) +
  geom_point(shape = 16, size = 4, alpha = 0.6) +
  #facet_wrap(~ sampling.period, nrow = 2) + 
  scale_x_continuous(breaks = c(22, 24, 26, 28, 30, 32)) + 
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold")) + 
  labs(x = expression(paste("Temperature ",degree,"C")), y = expression(z-transformed~Metabolic~rate~(VCO[2]~min^{-1})))
dev.off()

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

series.mass <- select(dat, id, samp_period, orig_lizmass, inverseK_prior_temp2) %>% 
  group_by(id, samp_period) %>% 
  summarise(mean_mass = mean(orig_lizmass),
            mean_prior_temp = mean(inverseK_prior_temp2))

series.mass$z.log.mass <- scale(log(series.mass$mean_mass))
names(series.mass) <- c("Lizard", "sampling.period", "mean_mass", "z.log.mass", "prior_temp")
names(cor.int.slop)

cor.int.slop2 <- left_join(cor.int.slop, series.mass)
str(cor.int.slop2)

slopes <- filter(cor.int.slop, sampling.period == '5') %>% select(Ind.slope, Series.slope) %>% mutate(Slope = Ind.slope + Series.slope) %>% data.frame()
  
ggplot(data = slopes, aes(Slope_rounded)) + 
  stat_density(geom = "line") + 
  stat_bin(bins = 20, binwidth = 0.03, alpha = 0.6) +
  theme_bw() +
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold")) + 
  labs(x = "Individual slopes", y = "Frequency")

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


#relationship with mass
cor.int.slop2 %>% 
  ggplot(aes(y = Series.int, x = Series.slope, color = z.log.mass)) + 
  geom_errorbar(aes(x = Series.slope, ymin = lower.Series.int, ymax = upper.Series.int), width = 0, size = 0.8) + 
  geom_errorbarh(aes(x = Series.slope, xmin = lower.Series.slope, xmax = upper.Series.slope), size = 0.8) + 
  geom_point(shape = 1, fill = "white", size = 1, color = "black")  + 
  facet_wrap(~ Lizard) +  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = "Slope for Series", y = "Intercept for Series") +
  scale_colour_gradient(low = "paleturquoise", high = "paleturquoise4", guide = "colourbar")
#dev.off()

#Relationshp with prior temp
cor.int.slop2 %>% 
  ggplot(aes(y = Series.int, x = Series.slope, color = prior_temp)) + 
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
  m4 <- readRDS("output/rds/m4_sp_log")
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
plot(m4.V1)
summary(m4.V1)
posterior.mode(m4.VCV)
HPDinterval(as.mcmc(m4.VCV))

#VCV matrix for Between ID - uses cor_cov matrix
# B = posterior.mode

mat <-posterior.mode(m4.VCV)
loctions3 <- grep("id", names(mat))
B <- mat[loctions3]

cor_cov_matrices(B = B, names = c(sort(unique(data$incb_temp))))

#Between individual correlations - posterior modes
write.csv(cor_cov_matrices(B = B, names = c(sort(unique(dat$incb_temp))))$cor, "output/table/betweenID_correlations.csv")
write.csv(cor_cov_matrices(B = B, names = c(sort(unique(dat$incb_temp))))$cov, "output/table/betweenID_covariance.csv")

#Between individual correlations - lower credible intervals
lower_mat <- HPDinterval(as.mcmc(m4.VCV))[,1]
locs_bw_lower <- grep("id", names(lower_mat))
low_B <- lower_mat[locs_bw_lower]

cor_cov_matrices(B = low_B, names = c(sort(unique(dat$incb_temp))))
write.csv(cor_cov_matrices(B = low_B, names = c(sort(unique(dat$incb_temp))))$cov,"output/table/betweenID_covariance_lower.csv")
write.csv(cor_cov_matrices(B = low_B, names = c(sort(unique(dat$incb_temp))))$cor,"output/table/betweenID_correlation_lower.csv")

#Between individual correlations - upper credible intervals
upper_mat <- HPDinterval(as.mcmc(m4.VCV))[,2]
locs_bw_upper <- grep("id", names(upper_mat))
upper_B <- upper_mat[locs_bw_upper]

cor_cov_matrices(B = upper_B, names = c(sort(unique(dat$incb_temp))))
write.csv(cor_cov_matrices(B = upper_B, names = c(sort(unique(dat$incb_temp))))$cov,"output/table/betweenID_covariance_upper.csv")
write.csv(cor_cov_matrices(B = upper_B, names = c(sort(unique(dat$incb_temp))))$cor,"output/table/betweenID_correlation_upper.csv")

#Within individual correlations
mat <-posterior.mode(m4.VCV)
loctions_w <- grep("units", names(mat))
w <- mat[loctions_w]

cor_cov_matrices(B = w, names = c(sort(unique(dat$incb_temp))))
write.csv(cor_cov_matrices(B = w, names = c(sort(unique(dat$incb_temp))))$cor, "output/table/withinID_correlations.csv")
write.csv(cor_cov_matrices(B = w, names = c(sort(unique(dat$incb_temp))))$cov, "output/table/withinID_covariance.csv")

#Within individual correlations - lower credible intervals
locs_wi_lower <- grep("units", names(lower_mat))
low_W <- lower_mat[locs_wi_lower]

cor_cov_matrices(B = low_W, names = c(sort(unique(dat$incb_temp))))
write.csv(cor_cov_matrices(B = low_W, names = c(sort(unique(dat$incb_temp))))$cor, "output/table/withinID_lower_correlations.csv")
write.csv(cor_cov_matrices(B = low_W, names = c(sort(unique(dat$incb_temp))))$cov, "output/table/withinID_lower_covariance.csv")

#Within individual correlations - upper credible intervals
locs_w_upper <- grep("units", names(upper_mat))
upper_W <- upper_mat[locs_w_upper]

cor_cov_matrices(B = upper_W, names = c(sort(unique(dat$incb_temp))))
write.csv(cor_cov_matrices(B = upper_W, names = c(sort(unique(dat$incb_temp))))$cor, "output/table/withinID_upper_correlations.csv")
write.csv(cor_cov_matrices(B = upper_W, names = c(sort(unique(dat$incb_temp))))$cov, "output/table/withinID_upper_covariance.csv")

#################################################################
#Graphing BETWEEN individual correlation matrix (posterior mode)#
#################################################################
#Because the matrix is not postive definitive... lets hand calculate this shit 
#Cor(x,y) = cov(x,y)/(sd(x) * sd(y))
#Need to use covariance matrix 

between.cov <- cor_cov_matrices(B = B, names = c(sort(unique(dat$incb_temp))))$cov
between.cor <- matrix(nrow = 6, ncol = 6)
row.names(between.cor) <- x
colnames(between.cor) <- x

#diagonals = 1 
diag(between.cor) <- rep(1, 6)

#22, 24
between.cor[2, 1] <- between.cov[2,1] / (sqrt(between.cov[1,1])*sqrt(between.cov[2,2]))

#22, 26
between.cor[3, 1] <- between.cov[3,1] / (sqrt(between.cov[1,1])*sqrt(between.cov[3,3]))

#22, 28
between.cor[4, 1] <- between.cov[4,1] / (sqrt(between.cov[1,1])*sqrt(between.cov[4,4]))

#22, 30
between.cor[5, 1] <- between.cov[5,1] / (sqrt(between.cov[1,1])*sqrt(between.cov[5,5]))

#22, 32
between.cor[6, 1] <- between.cov[6,1] / (sqrt(between.cov[1,1])*sqrt(between.cov[6,6]))

#posterior mode - upper triangle only
bw.cor <- cor_cov_matrices(B = B, names = c(sort(unique(dat$incb_temp))))$cor
melted.between_cor <- melt(get_upper_tri(bw.cor))
names(melted.between_cor)[3] <- "Correlation"
melted.between_cor$Correlation[c(1, 8, 15, 22, 29, 36)] <- rep(c(1),6)

#lower
bw.cor_lower <- cor_cov_matrices(B = low_B, names = c(sort(unique(dat$incb_temp))))$cor
melted.bw.cor.lower <- melt(get_lower_tri(bw.cor_lower))
names(melted.bw.cor.lower)[3] <- "Lower"
melted.between_cor$Lower <- melted.bw.cor.lower$Lower

#upper
bw.cor_upper <- cor_cov_matrices(B = upper_B, names = c(sort(unique(dat$incb_temp))))$cor
melted.bw.cor.upper <- melt(get_lower_tri(bw.cor_upper))
names(melted.bw.cor.upper)[3] <- "Upper"
melted.between_cor$Upper <- melted.bw.cor.upper$Upper

melted.between_cor$Correlation_labs <- paste0(melted.between_cor$Correlation, "\n" , " (", melted.between_cor$Lower, ",", melted.between_cor$Upper, ") ")

melted.between_cor$Correlation_labs[c(7, 13, 14, 19, 20, 21, 25, 26, 27, 28, 31, 32, 33, 34, 35)] <- c(" ")
melted.between_cor$Correlation_labs[c(1, 8, 15, 22, 29, 36)] <- c("1")

#PLOT FOR BETWEEN ID correlation
#pdf("output/fig/betweenID_cor.pdf", 7, 7)
ggplot(data = melted.between_cor, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile(colour = "White") + 
  scale_fill_continuous(low = "white", high = "red2",
                        limits = c(-1,1),
                        labels = c(-1, -0.5, 0, 0.5, 1),
                        breaks = c(-1, -0.5, 0, 0.5, 1), na.value = "white") +
  scale_x_continuous(breaks = c(sort(unique(dat$incb_temp)))) + 
  scale_y_continuous(breaks = c(sort(unique(dat$incb_temp)))) +
  #geom_text(aes(Var1, Var2, label = Correlation_labs), color = "white", size = 3) +
  labs(title = "Between ID correlations",
       subtitle = "Shrinkage method used") + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = "none",
        legend.direction = "horizontal") 
  #guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) 

#dev.off()

#################################################################
#Graphing BETWEEN individual covariance matrix (posterior mode)#
#################################################################

#posterior mode
bw.cov <- cor_cov_matrices(B = B, names = c(sort(unique(dat$incb_temp))))$cov
get_lower_tri(bw.cov)
melted.between_cov <- melt(get_lower_tri(bw.cov))
names(melted.between_cov)[3] <- "Covariance"

#lower
bw.cov_lower <- cor_cov_matrices(B = low_B, names = c(sort(unique(dat$incb_temp))))$cov
get_lower_tri(bw.cov_lower)
melted.bw.cov.lower <- melt(get_lower_tri(bw.cov_lower))
names(melted.bw.cov.lower)[3] <- "Lower"
melted.between_cov$Lower <- melted.bw.cov.lower$Lower

#upper
bw.cov_upper <- cor_cov_matrices(B = upper_B, names = c(sort(unique(dat$incb_temp))))$cov
melted.bw.cov.upper <- melt(get_lower_tri(bw.cov_upper))
names(melted.bw.cov.upper)[3] <- "Upper"
melted.between_cov$Upper <- melted.bw.cov.upper$Upper

melted.between_cov$Covariance_lab <- paste0(melted.between_cov$Covariance, "\n" , " (", melted.between_cov$Lower, ",", melted.between_cov$Upper, ") ")
melted.between_cov$Covariance_lab[c(7, 13, 14, 19, 20, 21, 25, 26, 27, 28, 31, 32, 33, 34, 35)] <- " "

#PLOT FOR BETWEEN ID covariance
#pdf("output/fig/betweenID_cov.pdf", 7, 7)
ggplot(data = melted.between_cov, aes(x = Var1, y = Var2, fill = Covariance)) +
  geom_tile(colour = "white") + 
  scale_fill_continuous(low = "navyblue", high = "orangered2", na.value = 'White') +
  scale_x_continuous(breaks = c(sort(unique(dat$incb_temp)))) + 
  scale_y_continuous(breaks = c(sort(unique(dat$incb_temp)))) +
  geom_text(aes(Var1, Var2, label = Covariance_lab), color = "white", size = 3) +
  labs(title = "Between ID covariance") + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.3,0.85),
        legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                title.position = "top", title.hjust = 0.5))
#dev.off()

#################################################################
#Graphing WITHIN individual correlation matrix (posterior mode)##
#################################################################
#posterior mode
w.cor <- cor_cov_matrices(B = w, names = c(sort(unique(dat$incb_temp))))$cor
melted.within_cor <- melt(get_lower_tri(w.cor))
names(melted.within_cor)[3] <- "Correlation"

#lower
w.cor_lower <- cor_cov_matrices(B = low_W, names = c(sort(unique(dat$incb_temp))))$cor
melted.w.cor.lower <- melt(get_lower_tri(w.cor_lower))
names(melted.w.cor.lower)[3] <- "Lower"
melted.within_cor$Lower <- melted.w.cor.lower$Lower

#upper
w.cor_upper <- cor_cov_matrices(B = upper_W, names = c(sort(unique(dat$incb_temp))))$cor
melted.w.cor.upper <- melt(get_lower_tri(w.cor_upper))
names(melted.w.cor.upper)[3] <- "Upper"
melted.within_cor$Upper <- melted.w.cor.upper$Upper

melted.within_cor$Correlation_labs <- paste0(melted.within_cor$Correlation, "\n" , " (", melted.within_cor$Lower, ",", melted.within_cor$Upper, ") ")
melted.within_cor$Correlation_labs[c(1, 8, 15, 22, 29, 36)] <- rep(c(1),6)
melted.within_cor$Correlation_labs[c(7, 13, 14, 19, 20, 21, 25, 26, 27, 28, 31, 32, 33, 34, 35)] <- " "

#PLOT FOR WITHIN ID correlation
#pdf("output/fig/withinID_cor.pdf", 7, 7)
ggplot(data = melted.within_cor, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile(colour = "white") + 
  scale_fill_continuous(low = "navyblue", high = "orangered2",
                        limits = c(-1,1),
                        labels = c(-1, -0.5, 0, 0.5, 1),
                        breaks = c(-1, -0.5, 0, 0.5, 1), na.value = "white") +
  scale_x_continuous(breaks = c(sort(unique(dat$incb_temp)))) + 
  scale_y_continuous(breaks = c(sort(unique(dat$incb_temp)))) +
  geom_text(aes(Var1, Var2, label = Correlation_2), color = "white", size = 3) +
  labs(title = "Within ID correlations") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.4,0.85),
        legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

#dev.off()

#################################################################
#Graphing WITHIN individual covariance matrix (posterior mode)###
#################################################################

#posterior mode
w.cov <- cor_cov_matrices(B = w, names = c(sort(unique(dat$incb_temp))))$cov
melted.within_cov <- melt(get_lower_tri(w.cov))
names(melted.within_cov)[3] <- "Covariance"

#lower
w.cov_lower <- cor_cov_matrices(B = low_W, names = c(sort(unique(dat$incb_temp))))$cov
melted.w.cov.lower <- melt(get_lower_tri(w.cov_lower))
names(melted.w.cov.lower)[3] <- "Lower"
melted.within_cov$Lower <- melted.w.cov.lower$Lower

#upper
w.cov_upper <- cor_cov_matrices(B = upper_W, names = c(sort(unique(dat$incb_temp))))$cov
melted.w.cov.upper <- melt(get_lower_tri(w.cov_upper))
names(melted.w.cov.upper)[3] <- "Upper"
melted.within_cov$Upper <- melted.w.cov.upper$Upper

melted.within_cov$Covariance_labs <- paste0(melted.within_cov$Covariance, "\n" , " (", melted.within_cov$Lower, ",", melted.within_cov$Upper, ") ")
melted.within_cov$Covariance_labs[c(7, 13, 14, 19, 20, 21, 25, 26, 27, 28, 31, 32, 33, 34, 35)] <- " "


#PLOT FOR WITHIN ID covariance
#pdf("output/fig/withinID_cov.pdf", 7, 7)
ggplot(data = melted.within_cov, aes(x = Var1, y = Var2, fill = Covariance)) +
  geom_tile(colour = "white") + 
  scale_fill_continuous(low = "navyblue", high = "orangered2", na.value = "white") +
  scale_x_continuous(breaks = c(sort(unique(dat$incb_temp)))) + 
  scale_y_continuous(breaks = c(sort(unique(dat$incb_temp)))) +
  geom_text(aes(Var1, Var2, label = Covariance_labs), color = "white", size = 3) +
  labs(title = "Within ID covariance") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.4,0.85),
        legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
#dev.off()

