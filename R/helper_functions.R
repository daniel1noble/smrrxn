#function to unpack chains into usable VCV

VCV.unpack <- function(dir.rds = "~/Dropbox/smrrxn/output/rds/m1.log_t22"){
  my.output <- readRDS(dir.rds)
  my.VCV <- lapply(my.output, function(m) m$VCV)
  my.chains <- do.call(rbind, my.VCV)
  return(my.chains)
}

#function to unpack chains into usable VCV
Sol.unpack <- function(dir.rds = "~/Dropbox/smrrxn/output/rds/m1.log_t22"){
  my.output <- readRDS(dir.rds)
  my.Sol <- lapply(my.output, function(m) m$Sol)
  my.chains <- do.call(rbind, my.Sol)
  return(my.chains)
}

#function to run diagnostics of VCV chains
VCV.check <- function(dir.rds = "~/Dropbox/smrrxn/output/rds/m1.log_t22"){
  my.output <- readRDS(dir.rds)
  my.VCV <- lapply(my.output, function(m) m$VCV)
  my.chains <- do.call(rbind, my.VCV)
  my.3chains <- do.call(coda::mcmc.list, my.VCV)
  plot(my.3chains)
  coda::autocorr.plot(as.mcmc(my.chains))
  coda::gelman.diag(my.3chains, multivariate = F)
}

#function to run diagnostics of VCV chains
Sol.check <- function(dir.rds = "~/Dropbox/smrrxn/output/rds/m1.log_t22"){
  my.output <- readRDS(dir.rds)
  my.Sol <- lapply(my.output, function(m) m$Sol)
  my.chains <- do.call(rbind, my.Sol)
  my.3chains <- do.call(coda::mcmc.list, my.Sol)
  plot(my.3chains)
  coda::autocorr.plot(as.mcmc(my.chains))
  coda::gelman.diag(my.3chains, multivariate = F)
}

my.DIC <-function(dir.rds = "~/Dropbox/smrrxn/output/rds/m1.log"){
  my.output <- readRDS(dir.rds)
  lapply(my.output, function(m) m$DIC)
}


rpt.Slope.mode <- function(VCVpooledchains, idslope, seriesslope){
  r.slope <- VCVpooledchains[,idslope] / ( VCVpooledchains[,idslope] + VCVpooledchains[,seriesslope]) 
  
  slope.tab <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(slope.tab) <- c("estimate","lower", "upper")
  rownames(slope.tab) <- c("rpt.slope")
  
  slope.tab[1,1] <- MCMCglmm::posterior.mode(r.slope)
  slope.tab[1,2:3] <- coda::HPDinterval(as.mcmc(r.slope))
  
  return(slope.tab)
}

rpt.Slope.mean <- function(VCVpooledchains, idslope, seriesslope){
  r.slope <- VCVpooledchains[,idslope] / ( VCVpooledchains[,idslope] + VCVpooledchains[,seriesslope]) 
  
  slope.tab <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(slope.tab) <- c("estimate","lower", "upper")
  rownames(slope.tab) <- c("rpt.slope")
  
  slope.tab[1,1] <- mean(r.slope)
  slope.tab[1,2:3] <- coda::HPDinterval(as.mcmc(r.slope))
  
  return(slope.tab)
}

FV.rpt.Int.mode <- function(VCVpooledchains){
  id.intercept <- colnames(VCVpooledchains)[1]
  series.intercept <- colnames(VCVpooledchains)[5]
  
  r.int <- VCVpooledchains[,id.intercept] / (VCVpooledchains[,id.intercept] + VCVpooledchains[,series.intercept])
  
  int.tab <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(int.tab) <- c("estimate","lower", "upper")
  rownames(int.tab) <- c("rpt.int")
  
  int.tab[1,1] <- MCMCglmm::posterior.mode(r.int)
  int.tab[1,2:3] <- coda::HPDinterval(as.mcmc(r.int))
  
  return(int.tab)
}

FV.rpt.Int.mean <- function(VCVpooledchains){
  id.intercept <- colnames(VCVpooledchains)[1]
  series.intercept <- colnames(VCVpooledchains)[5]
  
  r.int <- VCVpooledchains[,id.intercept] / (VCVpooledchains[,id.intercept] + VCVpooledchains[,series.intercept])
  
  int.tab <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(int.tab) <- c("estimate","lower", "upper")
  rownames(int.tab) <- c("rpt.int")
  
  int.tab[1,1] <- mean(r.int)
  int.tab[1,2:3] <- coda::HPDinterval(as.mcmc(r.int))
  
  return(int.tab)
}

FV.rpt.Int.long.mode <- function(VCVpooledchains){
  id.intercept <- colnames(VCVpooledchains)[1]
  series.intercept <- colnames(VCVpooledchains)[5]
  
  r.int.long <- VCVpooledchains[,id.intercept] / (VCVpooledchains[,id.intercept] + VCVpooledchains[,series.intercept] + VCVpooledchains[,"units"])
  
  int.tab <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(int.tab) <- c("estimate","lower", "upper")
  rownames(int.tab) <- c("rpt.int.long")
  
  int.tab[1,1] <- MCMCglmm::posterior.mode(r.int.long)
  int.tab[1,2:3] <- coda::HPDinterval(as.mcmc(r.int.long))
  
  return(int.tab)
}

FV.rpt.Int.long.mean <- function(VCVpooledchains){
  id.intercept <- colnames(VCVpooledchains)[1]
  series.intercept <- colnames(VCVpooledchains)[5]
  
  r.int.long <- VCVpooledchains[,id.intercept] / (VCVpooledchains[,id.intercept] + VCVpooledchains[,series.intercept] + VCVpooledchains[,"units"])
  
  int.tab <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(int.tab) <- c("estimate","lower", "upper")
  rownames(int.tab) <- c("rpt.int.long")
  
  int.tab[1,1] <- mean(r.int.long)
  int.tab[1,2:3] <- coda::HPDinterval(as.mcmc(r.int.long))
  
  return(int.tab)
}

FV.rpt.Int.short.mode <- function(VCVpooledchains){
  id.intercept <- colnames(VCVpooledchains)[1]
  series.intercept <- colnames(VCVpooledchains)[5]
  
  r.int.short <- (VCVpooledchains[,id.intercept] + VCVpooledchains[,series.intercept]) / (VCVpooledchains[,id.intercept] + VCVpooledchains[,series.intercept] + VCVpooledchains[,"units"])
  
  int.tab <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(int.tab) <- c("estimate","lower", "upper")
  rownames(int.tab) <- c("rpt.int.short")
  
  int.tab[1,1] <- MCMCglmm::posterior.mode(r.int.short)
  int.tab[1,2:3] <- coda::HPDinterval(as.mcmc(r.int.short))
  
  return(int.tab)
}

FV.rpt.Int.short.mean <- function(VCVpooledchains){
  id.intercept <- colnames(VCVpooledchains)[1]
  series.intercept <- colnames(VCVpooledchains)[5]
  
  r.int.short <- (VCVpooledchains[,id.intercept] + VCVpooledchains[,series.intercept]) / (VCVpooledchains[,id.intercept] + VCVpooledchains[,series.intercept] + VCVpooledchains[,"units"])
  
  int.tab <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(int.tab) <- c("estimate","lower", "upper")
  rownames(int.tab) <- c("rpt.int.short")
  
  int.tab[1,1] <- mean(r.int.short)
  int.tab[1,2:3] <- coda::HPDinterval(as.mcmc(r.int.short))
  
  return(int.tab)
}

SW.rpt.Temp.mode <- function(VCVpooledchains, temp = log(22)){
  id.intercept <- colnames(VCVpooledchains)[1]
  id.slope <- colnames(VCVpooledchains)[4]
  int.slope.covar <- colnames(VCVpooledchains)[2]
  temp <- temp
  
  rpt.temp <- (VCVpooledchains[,id.intercept] + (2 * VCVpooledchains[,int.slope.covar] * temp) + (VCVpooledchains[,id.slope] * temp^2)) / (VCVpooledchains[,id.intercept] + (2 * VCVpooledchains[,int.slope.covar] * temp) + (VCVpooledchains[,id.slope] * temp^2) + VCVpooledchains[,"units"])
  
  rpt.tab <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(rpt.tab) <- c("estimate","lower", "upper")
  rownames(rpt.tab) <- c(paste0("rpt", ".", exp(temp)))
  
  rpt.tab[1,1] <- MCMCglmm::posterior.mode(rpt.temp)
  rpt.tab[1,2:3] <- coda::HPDinterval(as.mcmc(rpt.temp))
  
  return(rpt.tab)
} 

SW.rpt.Temp.mean <- function(VCVpooledchains, temp = log(22)){
  id.intercept <- colnames(VCVpooledchains)[1]
  id.slope <- colnames(VCVpooledchains)[4]
  int.slope.covar <- colnames(VCVpooledchains)[2]
  temp <- temp
  
  rpt.temp <- (VCVpooledchains[,id.intercept] + (2 * VCVpooledchains[,int.slope.covar] * temp) + (VCVpooledchains[,id.slope] * temp^2)) / (VCVpooledchains[,id.intercept] + (2 * VCVpooledchains[,int.slope.covar] * temp) + (VCVpooledchains[,id.slope] * temp^2) + VCVpooledchains[,"units"])
  
  rpt.tab <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(rpt.tab) <- c("estimate","lower", "upper")
  rownames(rpt.tab) <- c(paste0("rpt", ".", exp(temp)))
  
  rpt.tab[1,1] <- mean(rpt.temp)
  rpt.tab[1,2:3] <- coda::HPDinterval(as.mcmc(rpt.temp))
  
  return(rpt.tab)
} 

CS.rpt.Temp <- function(VCVpooledchains, temp = 22){
  temp <- c(temp)
  var.temp <- paste0("t_",temp)
  btw_var <- paste0("trait",var.temp,":","trait",var.temp,".id")
  resid_var <- paste0("trait",var.temp,":","trait",var.temp,".units")
  
  cs.rpt.temp <- VCVpooledchains[,btw_var] / ( VCVpooledchains[,btw_var] + VCVpooledchains[,resid_var])
  
  rpt.tab <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(rpt.tab) <- c("estimate","lower", "upper")
  rownames(rpt.tab) <- c(paste0("rpt", ".", temp))
  
  rpt.tab[1,1] <- MCMCglmm::posterior.mode(cs.rpt.temp)
  rpt.tab[1,2:3] <- coda::HPDinterval(as.mcmc(cs.rpt.temp))
  
  return(rpt.tab)
}


my.cor.cov.matrices <- function(VCVpooledchains, type = "between"){
  matrix <- posterior.mode(VCVpooledchains)
  if(type == "between"){
  locations <- grep("id", names(matrix))
  } else{
    locations <- grep("units", names(matrix))
  }
  VCV.matrix.all <- matrix[locations]
  VCovM <- matrix(VCV.matrix.all, nrow = 6, ncol = 6)
  if(corpcor::is.positive.definite(VCovM)){
    VCorM <- cov2cor(VCovM)
  } else{
    VCorM <- cov2cor(corpcor::make.positive.definite(VCovM))
    warning("Covariance matrix not positive-definitive, estimates may not be accurate")
  }
  
  colnames(VCorM) <- colnames(VCovM) <- seq(22,32, by = 2)
  rownames(VCorM) <- rownames(VCovM) <- seq(22,32, by = 2)
  
  return(list(cov = VCovM, cor = VCorM))
}

my.postmode.HPD <- function(pooledchains){
  results <- cbind(posterior.mode(pooledchains), HPDinterval(as.mcmc(pooledchains)))
  colnames(results)[1] <- "mode"
  return(results)
}

my.postmean.HPD <- function(pooledchains){
  results <- cbind(colMeans(pooledchains), HPDinterval(as.mcmc(pooledchains)))
  colnames(results)[1] <- "mean"
  return(results)
}

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
      (incubation.temperatures[i] * (post$inverseK_incb_temp + post[, names(post) == lizard.slope] + post[, names(post) == sampling.slope]))  # temp * (global.slope +  slope for THIS lizard in THIS sampling period)
    # we leave out mass, since average is zero, so we are calcualting this for an average log(mass) liard
  }))
  
  #compiling predictions into data frame
  data.frame(predicted = c(out),
             Temperature = rep(incubation.temperatures, each = nrow(out)),
             Lizard = id, 
             stringsAsFactors = FALSE)
}

#Tabulationg FV temp models
tab_FV_temp.mode <- function(pooledchainsSol, pooledchainVCV){
  #Setting up table for all temp-centered models
  tab <- as.data.frame(matrix(ncol = 3, nrow = length(colnames(pooledchainsSol)) + 7))
  colnames(tab) <- c("estimate", "lower", "upper")
  rownames(tab) <- c(names(posterior.mode(pooledchainsSol)), names(posterior.mode(pooledchainVCV))[c(1,2,4,5,7,8,9)])
  
  tab[1:4,] <- my.postmode.HPD(pooledchainsSol)
  tab[5:11,] <- my.postmode.HPD(pooledchainVCV)[c(1,2,4,5,7,8,9),]
  
  return(tab)
}

tab_FV_temp.mean <- function(pooledchainsSol, pooledchainVCV){
  #Setting up table for all temp-centered models
  tab <- as.data.frame(matrix(ncol = 3, nrow = length(colnames(pooledchainsSol)) + 7))
  colnames(tab) <- c("estimate", "lower", "upper")
  rownames(tab) <- c(names(posterior.mode(pooledchainsSol)), names(posterior.mode(pooledchainVCV))[c(1,2,4,5,7,8,9)])
  
  tab[1:4,] <- my.postmean.HPD(pooledchainsSol)
  tab[5:11,] <- my.postmean.HPD(pooledchainVCV)[c(1,2,4,5,7,8,9),]
  
  return(tab)
}


#Tabulationg tab_reSW_FV_temp temp models
tab_SW_FV_temp.mode <- function(pooledchainsSol, pooledchainVCV){
  #Setting up table for all temp-centered models
  tab <- as.data.frame(matrix(ncol = 3, nrow = 9))
  colnames(tab) <- c("estimate", "lower", "upper")
  rownames(tab) <- c(names(posterior.mode(pooledchainsSol)), names(posterior.mode(pooledchainVCV))[c(1,3,4,5)])
  
  tab[1:5,] <- my.postmode.HPD(pooledchainsSol)
  tab[6:9,] <- my.postmode.HPD(pooledchainVCV)[c(1,3,4,5),]
  
  return(tab)
}

tab_SW_FV_temp.mean <- function(pooledchainsSol, pooledchainVCV){
  #Setting up table for all temp-centered models
  tab <- as.data.frame(matrix(ncol = 3, nrow = 9))
  colnames(tab) <- c("estimate", "lower", "upper")
  rownames(tab) <- c(colnames(pooledchainsSol), colnames(pooledchainVCV)[c(1,3,4,5)])
  
  tab[1:5,] <- my.postmean.HPD(pooledchainsSol)
  tab[6:9,] <- my.postmean.HPD(pooledchainVCV)[c(1,3,4,5),]
  
  return(tab)
}


#Tabulationg CS_temp models
tab_CS_temp <- function(pooledchainsSol, pooledchainVCV){
  #Setting up table for all temp-centered models
  tab <- as.data.frame(matrix(ncol = 3, 
                              nrow = length(c(names(posterior.mode(pooledchainsSol)), names(posterior.mode(pooledchainVCV))))))
  colnames(tab) <- c("estimate", "lower", "upper")
  rownames(tab) <- c(names(posterior.mode(pooledchainsSol)), names(posterior.mode(pooledchainVCV)))
  
  tab[1:3,] <- my.postmode.HPD(pooledchainsSol)
  tab[4:nrow(tab),] <- my.postmode.HPD(pooledchainVCV)
  
  return(tab)
}


brms_rpt <- function(model.rds.name = brms.m4.usall, temp = 22){
  y <- posterior_samples(model.rds.name)
  temp <- temp
  
  Vbetween <- paste0("sd_id__t", temp, "_Intercept")
  Vsession <- paste0("sd_samp_period__t", temp, "_Intercept")
  Vresidual <- paste0("sigma_t", temp)
  R <- y[names(y) == Vbetween] /(y[names(y) == Vbetween] + y[names(y) == Vsession] + y[names(y) == Vresidual])
  
  rpt_tab <- as.data.frame(matrix(nrow = 1, ncol = 4))
  rownames(rpt_tab) <- paste0("R_t_",temp)
  colnames(rpt_tab) <- colnames(posterior_summary(as.matrix(R)))
  rpt_tab[1,1:4] <- posterior_summary(as.matrix(R))
  return(rpt_tab)
}
