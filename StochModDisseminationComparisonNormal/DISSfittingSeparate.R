
library(gridExtra)
library(ggplot2)
library(binom)
library(beepr)
library(parallel)
library(plyr)
library(here)
#****model code****
source(here("StochModDisseminationComparison//DISSmodelFunc.R"))
#***likelihood****
source(here("StochModDisseminationComparison//DISSfittingFuncsNLL.R"))

#******************data to fit to********************
competenceDat <- read.csv(here(".//Data//datCiotaOnyango.csv"))
competenceDat <- competenceDat[competenceDat$Ref %in% "Onyango 2020",]
bin <- binom.confint(x=competenceDat$NumInf,n=competenceDat$ITotal,method="exact")
competenceDat$meanInf <- bin$mean
competenceDat$lowerInf <- bin$lower
competenceDat$upperInf <- bin$upper
#****************************************************


#****allowing all parameters to vary****
start <- Sys.time()
  trace<-3
  
#log_muV 1.544027e-01 log_infRate1  2.338052e-09 log_infRate2  5.494484e-08  

  #********initial parameter values****
  init.pars.fit <- c(
    log_prodRate1=log(50)
    ,log_prodRate2=log(10)
    ,log_cellSpread1=log(10^-3.5)
    ,log_cellSpread2=log(10^-3.5)
    ,log_escapeRate1=log(0.05)
    ,log_escapeRate2=log(0.0001))  
  #********optimise*******
  optim.vals <- optim(par = init.pars.fit
                      , objFXN
                      , fixed.params = virus_params()
                      , dat = competenceDat
                      , control = list(trace = trace)
                      , method = "Nelder-Mead" # 
                      , hessian = T)
  
  end <- Sys.time()
  end-start
  # save fitted parameters
  #AIC
  -2*(-optim.vals$value) + 2*2
  
  saveRDS(optim.vals,"DissModelFitSepAllParms220728.rds")
  
  
  #test <- readRDS("INFModelFitSepParms220713.rds")
  MLEfits <- optim.vals$par 
  estMuV <- exp(MLEfits)[1]
  estinfRate <- exp(MLEfits)[2]

