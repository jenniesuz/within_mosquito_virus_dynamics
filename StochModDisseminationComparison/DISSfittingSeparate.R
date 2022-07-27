
library(gridExtra)
library(ggplot2)
library(binom)
library(beepr)
library(parallel)
library(plyr)
library(here)
#****model code****
source(here("StochModInfectionComparison//INFmodelFunc.R"))
#***likelihood****
source(here("StochModInfectionComparison//INFfittingFuncsNLL.R"))

#******************data to fit to********************
competenceDat <- read.csv(here(".//Data//datCiotaOnyango.csv"))
competenceDat <- competenceDat[competenceDat$Ref %in% "Onyango 2020",]
bin <- binom.confint(x=competenceDat$NumInf,n=competenceDat$ITotal,method="exact")
competenceDat$meanInf <- bin$mean
competenceDat$lowerInf <- bin$lower
competenceDat$upperInf <- bin$upper
#****************************************************


#************************parameters for two 'treatments'****************************
virus_params <- function(   muV = 0.1
                            ,infRate1 = 10^-7.5
                            ,infRate2 = 10^-7.5
                            ,prodRate1 = 1              # note these parameters don't matter for infection
                            ,prodRate2 = 1              # but will be used in the dissemination model
                            ,cellSpread1 = 10^-4        #
                            ,cellSpread2 = 10^-4        #
                            ,escapeRate1 = 0.05         #
                            ,escapeRate2 = 0.05         #
                            ,cMax = 400                 #
)
  return(as.list(environment()))
#*****************************************************************************
start <- Sys.time()
  trace<-3
  #********initial parameter values****
  init.pars.fit <- c(
    log_muV=log(0.1)
    ,log_infRate1=log(10^-7.5)
  )
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
  
  saveRDS(optim.vals,"INFModelFitSameParms220714.rds")
  
  
  #test <- readRDS("INFModelFitSepParms220713.rds")
  MLEfits <- optim.vals$par 
  estMuV <- exp(MLEfits)[1]
  estinfRate <- exp(MLEfits)[2]

