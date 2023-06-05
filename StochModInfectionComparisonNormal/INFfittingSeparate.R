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
competenceDat <- read.csv(here("StochModInfectionComparisonNormal//datCiotaOnyango.csv"))
competenceDat <- competenceDat[competenceDat$Ref %in% "Ciota 2017",]
bin <- binom.confint(x=competenceDat$NumInf,n=competenceDat$ITotal,method="exact")
competenceDat$meanInf <- bin$mean
competenceDat$lowerInf <- bin$lower
competenceDat$upperInf <- bin$upper
#****************************************************


# test for difference between species in dose-response curve using logistic regression
noSppMod <- glm(NumInf/ITotal ~ ConcMax,data=competenceDat,family="binomial",weights=ITotal)
SppMod <- glm(NumInf/ITotal ~ ConcMax+Moz,data=competenceDat,family="binomial",weights=ITotal)
SppModInt <- glm(NumInf/ITotal ~ ConcMax*Moz,data=competenceDat,family="binomial",weights=ITotal)
AIC(noSppMod)
AIC(SppMod)
AIC(SppModInt)



#**********************get approx starting parameters*********************

virus_params <- function(   muV = 0.1
                            ,infRate1 = 10^-8
                            ,prodRate1 = 10            # note these parameters don't matter for infection
                            ,cellSpread1 = 10^-4        
                            ,escapeRate1 = 0.05         #
                            ,cMax = 400                 #
)
  return(as.list(environment()))


testAeg <- repeatInfModel(x=virus_params(infRate1=10^-9)
                          ,startingVirus=10^competenceDat$ConcMax[competenceDat$Moz %in% "Ae. aegypti"]
)
modAeg <- glm(num/denom~conc,family="binomial",weights=denom,data=testAeg)
predAeg<-predict(modAeg,type="response",newdata=newData)
predAeg <- cbind.data.frame("conc"=seq(4,10,0.5),"predVals"=predAeg)
#*********** values for Aedes albopictus *******

#*
testAlb <- repeatInfModel(x=virus_params(infRate1=10^-7)
                          ,startingVirus=10^competenceDat$ConcMax[competenceDat$Moz %in% "Ae. albopictus"]
)
modAlb <- glm(num/denom~conc,family="binomial",weights=denom,data=testAlb)
predAlb<-predict(modAlb,type="response",newdata=newData)
predAlb <- cbind.data.frame("conc"=seq(4,10,0.5),"predVals"=predAlb)
#
ggplot(competenceDat) +
  geom_point(aes(x=ConcMax,y=meanInf,col=Moz)) +
  geom_line(data=predAeg,aes(x=conc,y=predVals)) +
  geom_line(data=predAlb,aes(x=conc,y=predVals))





#************************parameters for two 'treatments'****************************
virus_params <- function(   muV = 0.1
                            ,infRate1 = 10^-8
                            ,infRate2 = 10^-8
                            ,prodRate1 = 1              # note these parameters don't matter for infection
                            ,cellSpread1 = 10^-4        
                            ,escapeRate1 = 0.05         #
                            ,cMax = 400                 #
)
  return(as.list(environment()))
#*****************************************************************************
init.pars.fit <- c(
  log_muV=log(0.1)
  ,log_infRate1=log(10^-9)
  ,log_infRate2=log(10^-7)
)

objFXN(fit.params=init.pars.fit                                                          ## paramters to fit
                   , fixed.params =virus_params()                                      ## fixed paramters
                   , dat=competenceDat
                   , nSimulations=30)

#************curve******************
start <- Sys.time()
testLik30 <- lapply(10^-seq(7,9,0.1)
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

dat <- cbind.data.frame(beta=10^-seq(7,9,0.1),ll=testLik30)

ggplot(dat) +
  geom_line(aes(x=log10(beta),y=ll)) +
  ylim(0,25)



#******************************************

start <- Sys.time()
  trace<-3

  init.pars.fit <- c(
    log_infRate1=log(10^-9)
    ,log_infRate2=log(10^-7)
  )
  
  #********optimise*******
  optim.vals <- optim(par = init.pars.fit
                      , objFXN
                      , fixed.params = virus_params()
                      , dat = competenceDat
                      , nSimulations = 30
                      , control = list(trace = trace
                                       #,abstol=0.05
                                       #,reltol=0.05
                                       ,maxit=200
                                      )
                      , method ="SANN" #"Nelder-Mead" #"SANN" #
                     )
  
  end <- Sys.time()
  end-start

  -2*(-optim.vals$value) + 2*2

# epilson 0.05 SANN 1 run
 # $value
 # [1] 6.706336
#$par
# log_infRate1 log_infRate2 
# -20.50493    -16.78952   
  



  
  #test <- readRDS("INFModelFitSepParms220713.rds")
  MLEfits <- optim.vals$par 
  estMuV <- exp(MLEfits)[1]
  estinfRate <- exp(MLEfits)[2]

