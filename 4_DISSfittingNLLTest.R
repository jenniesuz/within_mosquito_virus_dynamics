library(gridExtra)
library(ggplot2)
library(binom)
library(beepr)

#********test fitting procedure and initial parameter value estimates

source("4_DISSfittingFuncsNLL.R")
source("4_DISSmodelFunc.R")
source("3_INFfittingNLLParameterEstimates.R")

#******************data to fit to********************
competenceDat <- read.csv("datCiotaOnyango.csv")
competenceDat <- competenceDat[competenceDat$Ref %in% "Onyango 2020",]
bin <- binom.confint(x=competenceDat$NumDiss,n=competenceDat$NumInf,method="exact")
competenceDat$meanDiss <- bin$mean
competenceDat$lowerDiss <- bin$lower
competenceDat$upperDiss <- bin$upper
#****************************************************


#*************************************************************
virus_params <- function(   muV=0.1
                            ,infRate=10^-7
                            ,prodRate =   10 
                            ,cellSpread = 10^-8
                            ,escapeRate = 0.005
                            ,cMax = 400            
                            ,hMax = 900             
)
  return(as.list(environment()))
#**************************************************************


testAlb <- sim.func(10^competenceDat$Conc.Min[1],parms=virus_params(infRate=albinfRate,muV=albMuV))

#**********first remove simulations where an infection didn't occur**********************
inf <- ddply(testAlb,.(run),summarise,occurred=sum(inf))  # establish if infection occurred
inf$occurred[inf$occurred>0] <- 1
noInf <- inf$run[inf$occurred == 0]                         # note run and concs where infection didn't occur

testAlb <- testAlb[!testAlb$run %in% noInf,]
#*******************************************************************************

testAlb$days <- round(testAlb$t/24,1)  # create a column of .1 days
diss2 <- dissSummaryFunc(testAlb)

ggplot(diss2) + 
  geom_point(data=competenceDat,aes(x=DPIDissorTrans,y=meanDiss))+
  geom_errorbar(data=competenceDat,aes(x=DPIDissorTrans,ymin=lowerDiss,ymax=upperDiss)) +
  geom_line(aes(x=time,y=proportionDisseminated)) +
  xlab("Time (days)") +
  ylab("Proportion of infected simulations with disseminated infection") +
  ylim(0,1) +
  labs(col=expression("Input virus concentration (log"[10]*")")) 


#***************************************************************************
init.pars.fit1 <- c(   
  log_prodRate =   log(10) 
  ,log_cellSpread = log(10^-8)  
  ,log_escapeRate = log(0.005)   
)
#**************Optimise*******************************************
trace <- 3

optim.vals <- optim(par = init.pars.fit1
                    , objFXN
                    , fixed.params = virus_params(infRate=albinfRate,muV=albMuV)
                    , dat = competenceDat
                    , control = list(trace = trace, maxit=100,reltol = 0.001)
                    , method = "Nelder-Mead" # 
                    , hessian = T)
optim.vals # convergence 0 means algorithm converged
beep()

MLEfits <- optim.vals$par 
exp(MLEfits)

##

test <- sim.func(competenceDat$Conc.Min[1],parms=virus_params(infRate=albinfRate,muV=albMuV
  ,prodRate=exp(MLEfits)[1]
  ,cellSpread=exp(MLEfits)[2]
  ,escapeRate=exp(MLEfits)[3]
))

tempDiss <- dissSummaryFunc(testAlb)


ggplot(tempDiss) + 
  geom_point(data=competenceDat,aes(x=DPIDissorTrans,y=meanDiss,fill=factor(Moz)),shape=21) +
  geom_errorbar(data=competenceDat,aes(x=DPIDissorTrans,ymin=lowerDiss,ymax=upperDiss)) +
  geom_line(aes(x=time,y=proportionDisseminated)) +
  xlab("Time (days)") +
  ylab("Proportion of infected simulations with disseminated infection") +
  ylim(0,1) +
  labs(col=expression("Input virus concentration (log"[10]*")")) 





