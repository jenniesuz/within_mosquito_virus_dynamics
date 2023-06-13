library(mvtnorm)
library(gridExtra)
library(ggplot2)
library(binom)
library(beepr)
library(parallel)
library(here)
#****model code****
source(here("StochModDisseminationComparisonNormal//InfDissModelFuncs.R"))
source(here("StochModDisseminationComparisonNormal//InfDissRepeatModelFunc.R"))
#***likelihood****
source(here("StochModDisseminationComparisonNormal//DISSfittingFuncsNLL.R"))

#******************data to fit to********************
competenceDat <- read.csv(here("StochModInfectionComparisonNormal//datCiotaOnyango.csv"))
competenceDat <- competenceDat[competenceDat$Ref %in% "Onyango 2020",]
bin <- binom.confint(x=competenceDat$NumDiss,n=competenceDat$NumInf,method="exact")
competenceDat$meanDiss <- bin$mean
competenceDat$lowerDiss <- bin$lower
competenceDat$upperDiss <- bin$upper
#****************************************************

#Fit GLM to data
modDatnoSpp <- glm(NumDiss/NumInf~DPIDissorTrans
              ,family="binomial"
              ,weights=NumInf
              ,data=competenceDat)
AIC(modDatnoSpp)
modDatSpp <- glm(NumDiss/NumInf~DPIDissorTrans + Moz
                   ,family="binomial"
                   ,weights=NumInf
                   ,data=competenceDat)
AIC(modDatSpp)

newData <- c("DPIDissorTrans"=c(1:14))
preds<-predict(modAeg,type="response",newdata=newData)
preds <- cbind.data.frame("DPIDissorTrans"=c(1:14),"predVals"=preds)

#******plot***
plotDat <- ggplot(competenceDat) +
  geom_point(aes(x=DPIDissorTrans,y=meanDiss,col=Moz)) +
  geom_errorbar(aes(x=DPIDissorTrans,ymin=lowerDiss,ymax=upperDiss,col=Moz)) +
  xlim(0,20)

plotDat

#************get ball park starting parameters***********
virus_params <- function(   muV = 0.1 #0.1
                            ,infRate = 10^-9 #10^-8.8
                            ,prodRate = 100#170         #
                            ,cellSpread = 0.005 # 0.00005       
                            ,escapeRate = 0.005 # 0.25#0.005  #
                            ,cMax = 400  
                            ,hMax = 900             #
)
  return(as.list(environment()))

test <- repeatModel(virus_params(),startingVirus=10^competenceDat$ConcMax[1])
test <- do.call(rbind,test)
# get rid of runs where midgut infection didn't occur
test <- test[test$inf %in% 1,]   # make sure calculating the probability of dissemination given infection
testNums <- dissSummaryFunc(test)
testNums$props <- testNums$numberRunsDisseminated/testNums$totalSize

plotDat +
  geom_line(data=testNums,aes(x=time,y=props)) 

#*******************parameters*****************
virus_params <- function(   muV = 0.1
                            ,infRate1 = 10^-9 # 2.338052e-09
                            ,infRate2 = 10^-7
                            ,prodRate1 = 170  
                            ,prodRate2 = 50
                            ,cellSpread1 = 0.00023 
                            ,cellSpread2 = 0.00045
                            ,escapeRate1 = 0.25 # 0.0005
                            ,escapeRate2 = 0.00024
                            ,cMax = 400  
                            ,hMax= 900
)
return(as.list(environment()))

nll.binom(parms=virus_params()
          ,dat=competenceDat
          ,nSims=30)








#****all three different between species****
start <- Sys.time()
  trace<-3
  
  #********initial parameter values****
  init.pars.fit <- c(
    log_prodRate1=log(170)  #170
    ,log_prodRate2=log(50)
    ,log_cellSpread1=log(0.00023)
    ,log_cellSpread2=log(0.00045)
    ,log_escapeRate1=log(0.25)
    ,log_escapeRate2=log(0.00024)
  )  
  #********optimise*******
  optim.vals <- optim(par = init.pars.fit
                      , objFXN
                      , fixed.params = virus_params(muV=exp(-3.089103),infRate1=exp(-21.729066),infRate2=exp(-17.346550))
                      , dat = competenceDat
                      , control = list(trace = trace, maxit=200)
                      , method = "SANN" # 
                      , hessian = T)
  
  end <- Sys.time()
  end-start
  # save fitted parameters
  #AIC
  -2*(-optim.vals$value) + 2*6
  
  
  
  
saveRDS(optim.vals,"DissModelFitSepAllParms230608.rds")
  
fit <- readRDS("DissModelFitSepAllParms230608.rds")
