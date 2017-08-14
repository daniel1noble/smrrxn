setwd("~/Dropbox/smrrxn/")

#clear envir
rm(list = ls())

#load library
library(MCMCglmm)
library(parallel)

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

#Rshort
R.short <- ( m1.VCV[,"(Intercept):(Intercept).id"] + m1.VCV[,"(Intercept):(Intercept).series"] ) / ( m1.VCV[,"(Intercept):(Intercept).id"] + m1.VCV[,"(Intercept):(Intercept).series"]  + m1.VCV[,"units"] )

Table1[14,1] <- posterior.mode(R.short) #R.short
Table1[14,2:3] <- HPDinterval(as.mcmc(R.short)) #R.short CIs

#Rlong
R.long <-  m1.VCV[,"(Intercept):(Intercept).id"] / ( m1.VCV[,"(Intercept):(Intercept).id"] + m1.VCV[,"(Intercept):(Intercept).series"]  + m1.VCV[,"units"] )

Table1[15,1] <-  posterior.mode(R.long) #Rlong
Table1[15,2:3] <- HPDinterval(as.mcmc(R.long)) #Rlong CIs

write.csv(Table1, "output/table/Table1.csv")
write.csv(round(Table1,2), "output/table/Table1_rounded.csv")

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

reaction.norms <- output %>% group_by(Temperature, Lizard, sampling.period) %>% summarise(posterior.mode = posterior.mode(as.mcmc(predicted)), lowerCI = as.numeric(HPDinterval(as.mcmc(predicted)))[1], upperCI = as.numeric(HPDinterval(as.mcmc(predicted)))[2])

str(reaction.norms)

#Plotting reactions
#reaction.norms %>% ggplot(aes(Lizard, posterior.mode)) + geom_hline(yintercept = 0, linetype = 2) + geom_errorbar(aes(ymin = lowerCI, ymax = upperCI), width = 0) + geom_point() + coord_flip() + facet_wrap(~Temperature) + ylab("Respiration rate")

#pdf("output/fig/reaction.norms.pdf", 10, 6)
reaction.norms %>% ggplot(aes(x = Temperature, y = posterior.mode, group = Lizard, color = Lizard)) + geom_line() + geom_point(shape = 1, fill = "white", size = 1, color = "black") + facet_wrap(~ sampling.period, nrow = 2) + theme_bw() + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "Temperature (1/K)", y = expression(Metabolic~rate~(CO[2]~min^{-1})))
#dev.off()

#pdf("output/fig/ID.rxn.norm.pdf", 9,9)
reaction.norms %>% ggplot(aes(x = Temperature, y = posterior.mode, group = sampling.period, color = sampling.period)) + geom_line() + geom_point(shape = 1, fill = "white", size = 1, color = "black") + facet_wrap(~ Lizard) + theme_bw() + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "Temperature (1/K)", y = expression(Metabolic~rate~(CO[2]~min^{-1})))
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

output <- do.call("rbind", lapply(1:10, function(i){
  do.call("rbind", lapply(lizard.names, get.intercept.slope, post=cbind(m2.Sol, m2.VCV), sampling.period = i)) %>% mutate(sampling.period = i) %>% arrange(sampling.period) %>% mutate(Lizard = factor(Lizard, levels = unique(Lizard))) 
}))

cor.int.slop <- output %>% group_by(Lizard, sampling.period) %>% summarise(global.int = posterior.mode(as.mcmc(global.int)),
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
                                                                           upper.Series.slope = HPDinterval(as.mcmc(series.slope))[,2])


#Plotting ints with slopes

#pdf("output/fig/covariance.ID.int.slope.pdf", 10, 6)
fig2a <- cor.int.slop %>% ggplot(aes(y = Ind.int, x = Ind.slope)) + geom_errorbar(aes(x = Ind.slope, ymin = lower.Ind.int, ymax = upper.Ind.int), color = "#999999", width = 0, size = 0.25)+ geom_errorbarh(aes(x = Ind.slope, xmin = lower.Ind.slope, xmax = upper.Ind.slope), color = "#999999", size = 0.25) + geom_point(shape = 1, fill = "white", size = 1, color = "black")  + facet_wrap(~ sampling.period, nrow = 2) + theme_bw() + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "Slope for ID", y = "Intercept for ID") 
#dev.off()

#pdf("output/fig/covariance.series.int.slope.pdf", 10, 6)
fig2b <- cor.int.slop %>% ggplot(aes(y = Series.int, x = Series.slope)) + geom_errorbar(aes(x = Series.slope, ymin = lower.Series.int, ymax = upper.Series.int), width = 0, size = 0.25, color = "#999999")+ geom_errorbarh(aes(x = Series.slope, xmin = lower.Series.slope, xmax = upper.Series.slope), size = 0.1, color = "#999999") + geom_point(shape = 1, fill = "white", size = 0.25, color = "black")  + facet_wrap(~ sampling.period, nrow = 2) +  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "Slope for Series", y = "Intercept for Series") 
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
cor.int.slop %>% ggplot(aes(y = Series.int, x = Series.slope, color = sampling.period)) + geom_errorbar(aes(x = Series.slope, ymin = lower.Series.int, ymax = upper.Series.int), width = 0, size = 0.8)+ geom_errorbarh(aes(x = Series.slope, xmin = lower.Series.slope, xmax = upper.Series.slope), size = 0.8) + geom_point(shape = 1, fill = "white", size = 1, color = "black")  + facet_wrap(~ Lizard) +  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "Slope for Series", y = "Intercept for Series") 
#dev.off()



#Intercept model for calculating PCV 
#m3
if(m3){
  m3 <- mclapply(1:3, function(i) {
    MCMCglmm(z.log.co2pmin ~ 1,
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
  m3 <- readRDS("output/rds/m3")
}

