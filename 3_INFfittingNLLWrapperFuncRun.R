library(gridExtra)
library(ggplot2)
library(binom)
library(beepr)
library(parallel)

source("3_INFfittingFuncsNLL.R")
source("3_INFmodelFunc.R")

virus_params <- function(   muV = 0.1
                            ,infRate = 10^-7.5
                            ,prodRate =   1    
                            ,cellSpread = 10^-4  
                            ,escapeRate = 0.05    
                            ,cMax = 400  
)
return(as.list(environment()))

#******************data to fit to********************
competenceDat <- read.csv("datCiotaOnyango.csv")
competenceDat <- competenceDat[competenceDat$Ref %in% "Ciota 2017",]
bin <- binom.confint(x=competenceDat$NumInf,n=competenceDat$ITotal,method="exact")
competenceDat$meanInf <- bin$mean
competenceDat$lowerInf <- bin$lower
competenceDat$upperInf <- bin$upper
#****************************************************


#**************************Function to wrap fitting and simulations for iteration*************
fitToDataFunc <- function(x
                          ,virusDecay=0.1
                          ,infRate=10^-9.5
                          ,dataSet=competenceDat[(competenceDat$Moz %in% "Ae. aegypti"),]){
  a<-x
  trace<-3
  #********initial parameter values****
  init.pars.fit <- c(
    log_muV=log(virusDecay)
    ,log_infRate=log(infRate)
  )
  #********optimise*******
  optim.vals <- optim(par = init.pars.fit
                      , objFXN
                      , fixed.params = virus_params()
                      , dat = dataSet
                      , control = list(trace = trace, maxit = 500, reltol = 10^-7.5)
                      , method = "Nelder-Mead" # 
                      , hessian = T)
  # save fitted parameters
  MLEfits <- optim.vals$par 
  estMuV <- exp(MLEfits)[1]
  estinfRate <- exp(MLEfits)[2]
  
  return(list(optim.vals,c(estMuV,estinfRate)))
}
#********************************************************************************


bothInffits <- mclapply(1,fitToDataFunc,dataSet=competenceDat,infRate=10^-8.5)

  
aegInffits <- mclapply(1:100,fitToDataFunc,dataSet=competenceDat[(competenceDat$Moz %in% "Ae. aegypti"),])
albInffits <- mclapply(1:100,fitToDataFunc,dataSet=competenceDat[(competenceDat$Moz %in% "Ae. albopictus"),],infRate=10^-7.5)

saveRDS(aegInffits,"aegInffits.rds")
saveRDS(albInffits,"albInffits.rds")








#*******************************Simulate with fitted values**********************
aeg <- readRDS("aegInffits.rds")
alb <- readRDS("albInffits.rds")

aegParams <- sapply(aeg,"[[",2)
albParams <- sapply(alb,"[[",2)

doses <- c(10^3,10^3.5,10^4,10^4.5,10^5,10^5.5,10^6,10^6.5,10^7,10^7.5,10^8,10^8.5,10^9,10^9.5,10^10)

#***************Ae aegypti*************
start_time <- Sys.time()

aegSims <- lapply(1:length(aegParams[2,]),function(x){
  sims <- mclapply(doses,
                         sim.func,parms=virus_params(muV=aegParams[1,x] ,infRate=aegParams[2,x]))
  doseSim <- do.call(rbind.data.frame,sims)
  infDat <- modOutFunc(doseSim)
  infDat$infRate <- aegParams[2,x]
  infDat$muV <- aegParams[1,x]
  return(infDat)
})

end_time <- Sys.time()
timeTaken <- end_time - start_time
timeTaken

saveRDS(aegSims,"aegInfSims.rds")

#*********Ae albopictus********

start_time <- Sys.time()

albSims <- lapply(1:length(albParams[2,]),function(x){
  sims <- mclapply(doses,
                   sim.func,parms=virus_params(muV=albParams[1,x] ,infRate=albParams[2,x]))
  doseSim <- do.call(rbind.data.frame,sims)
  infDat <- modOutFunc(doseSim)
  infDat$infRate <- albParams[2,x]
  infDat$muV <- albParams[1,x]
  return(infDat)
})

end_time <- Sys.time()
timeTaken <- end_time - start_time
timeTaken

saveRDS(albSims,"albInfSims.rds")

#********************************************************

