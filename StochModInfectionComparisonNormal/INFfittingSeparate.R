library(here)
library(binom)
library(ggplot2)
library(bbmle)
library(parallel)
library(adaptivetau)
library(mvtnorm)
#****model code****
source(here("StochModInfectionComparisonNormal//INFmodelFunc.R"))

source(here("StochModInfectionComparisonNormal//INFRepeatModelFunc.R"))
#***likelihood****
source(here("StochModInfectionComparisonNormal//INFfittingFuncsNLL.R"))

#******************data to fit to********************
competenceDat <- read.csv(here(".//Data//datCiotaOnyango.csv"))
competenceDat <- competenceDat[competenceDat$Ref %in% "Ciota 2017",]
bin <- binom.confint(x=competenceDat$NumInf,n=competenceDat$ITotal,method="exact")
competenceDat$meanInf <- bin$mean
competenceDat$lowerInf <- bin$lower
competenceDat$upperInf <- bin$upper
#****************************************************


#**********************get approx starting parameters*********************
ggplot(competenceDat) +
  geom_point(aes(x=ConcMax,y=meanInf,col=Moz))

virus_params <- function(   muV = 0.1
                            ,infRate1 = 10^-8
                            ,prodRate1 = 10            # note these parameters don't matter for infection
                            ,cellSpread1 = 10^-4        
                            ,escapeRate1 = 0.05         #
                            ,cMax = 400                 #
)
  return(as.list(environment()))

testAeg <- repeatInfModel(x=virus_params(infRate1=10^-8.5)
 ,startingVirus=10^competenceDat$ConcMax[competenceDat$Moz %in% "Ae. aegypti"]
)
testAeg

#************************parameters for two 'treatments'****************************
virus_params <- function(   muV = 0.1
                            ,infRate1 = 10^-8.5
                            ,infRate2 = 10^-8
                            ,prodRate1 = 10              # note these parameters don't matter for infection
                            ,cellSpread1 = 10^-4        
                            ,escapeRate1 = 0.05         #
                            ,cMax = 400                 #
)
  return(as.list(environment()))
#*****************************************************************************
init.pars.fit <- c(
  log_infRate1=log(10^-9)
  ,log_infRate2=log(10^-7.5)
)

objFXN(fit.params=init.pars.fit                                                          ## paramters to fit
                   , fixed.params =virus_params()                                      ## fixed paramters
                   , dat=competenceDat
                   , nSimulations=30)

#************curve******************
start <- Sys.time()
testLik30 <- lapply(10^-7#10^-seq(7,9,0.1)
                    ,function(x){
  loglik <- nll.binom(parms=virus_params(muV=0.1
                                         ,infRate1=x
                                         ,infRate2=x)
                      ,dat=competenceDat
                      ,nSims=30) 
  return(loglik)
})

testLik30 <- do.call(c,testLik30)
end <- Sys.time()
end-start

dat <- cbind.data.frame(beta=10^-seq(7,8,0.05),ll=testLik30)

ggplot(dat) +
  geom_line(aes(x=log10(beta),y=ll)) +
  ylim(0,25)




#******************************************

start <- Sys.time()
  trace<-3
  #********initial parameter values****
  init.pars.fit <- c(
    log_infRate1=log(10^-8)
    ,log_infRate2=log(10^-7.5)
  )
  #********optimise*******
  optim.vals <- optim(par = init.pars.fit
                      , objFXN
                      , fixed.params = virus_params()
                      , dat = competenceDat
                      , nSimulations = 30
                      , control = list(trace = trace
                                       ,abstol=0.001
                                       ,reltol=0.001)
                      , method = "Nelder-Mead" # 
                     )
  
  end <- Sys.time()
  end-start
  # save fitted parameters
  #AIC
  -2*(-optim.vals$value) + 2*2
  
  saveRDS(optim.vals,"INFModelFitSameParms230429.rds")
  
  
  #test <- readRDS("INFModelFitSepParms220713.rds")
  MLEfits <- optim.vals$par 
  estMuV <- exp(MLEfits)[1]
  estinfRate <- exp(MLEfits)[2]

