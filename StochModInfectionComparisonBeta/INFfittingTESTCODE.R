library(here)
library(binom)
library(bbmle)
library(emdbook)
library(EnvStats)
library(optimx)
library(ggplot2)
#****model code****
source(here("StochModInfectionComparisonBeta//INFRepeatModelFunc.R"))
#***likelihood****
source(here("StochModInfectionComparisonBeta//INFfittingFuncsNLLTEST.R"))

#************************parameters****************************
virus_params <- function(   muV = 0.1
                            ,infRate = 10^-7.5
                            ,prodRate =   1  
                            ,cellSpread = 10^-4  
                            ,escapeRate = 0.05
                            ,cMax = 400  
)
  return(as.list(environment()))
#*****************************************************************************

library(plyr)
#****generate data using 'true' parameter values****
# 30 mosquitoes, 3 virus concentrations
testDat <- repeatInfModel(x=virus_params()
                          ,startingVirus=c(10^5,10^6,10^6.5)
)
names(testDat) <- c("NumInf","ITotal","ConcMax")


#*****************************************************************************

testRun <- repeatInfModel(virus_params(infRate=10^-7.5
                                     ,muV=0.1)
                        ,startingVirus=c(10^4,10^5,10^6,10^7,10^8)
)

testRun$mean <- testRun$num/testRun$denom
names(testRun)[3]<- "ConcMax"


#***************explore likelihoods************
testLik30 <- nll.binom(parms=virus_params(muV=0.1
                              ,infRate=10^-7.5)
                              ,dat=testDat
                              ,nSims=30) 


testLik30 <- lapply(c(6,6.5,7,7.5),function(x){
  loglik <- nll.binom(parms=virus_params(muV=0.1
                               ,infRate=x)
            ,dat=testDat
            ,nSims=30) 
  return(loglik)
})

testLik30 <- do.call(c,testLik30)

# There is a trade-off between number of simulations used to
# estiamte the beta parameters and therefore accuracy of those
# estimations and the time it takes to fit the model


#********initial parameter values****
init.pars.fit <- c(
  log_infRate=log(10^-9)  
)
  
start <- Sys.time()

#********optimise*******
optim.vals <- optimx(par = init.pars.fit
                    , objFXN
                    , fixed.params = virus_params()
                    , nSimulations = 30
                    , dat = testDat
                    , control = list(trace = 3
                                     ,reltol=0.01)
                    , method = "Nelder-Mead" # 
                    , hessian = T)
  
 end <- Sys.time()
end-start

saveRDS(optim.vals,"OptimVals30simsInfRateOnly.RDS") # when 100 maxit takes c. 20-30 mins # when 200 sims took 1 hour



exp(coef(optim.vals))[1]
log10(exp(coef(optim.vals))[2])


#The Hessian matrix gives you the curvature of the likelihood function at the maximum likelihood
## estimate (MLE) of the fitted paramters. In other words, it tells you the second derivative around
## MLE, which can be used to estimate the covariance variance matrix of the MLE. This estimate of
## the covariance varance matrix is known as the Fisher information matrix and can be obtained by
## inverting the Hessian.
fisherInfMatrix <- solve(optim.vals$hessian) ## invert the Hessian, to estimate the covar-var matrix of parameter estimates





#*********plot
cint <- binom.confint(x=testDat$NumInf,n=testDat$ITotal,methods="exact")
testDat$mean <- cint$mean
testDat$lower <- cint$lower
testDat$upper <- cint$upper

plotDat <- ggplot(testDat) +
  geom_point(aes(x=ConcMax,y=mean)) +
  geom_errorbar(aes(x=ConcMax,ymin=lower,ymax=upper)) +
  ylim(0,1)
plotDat


fit30 <- repeatInfModel(virus_params(infRate=exp(coef(optim.vals))[2]
                                     ,muV=exp(coef(optim.vals))[1])
                        ,startingVirus=c(10^4,10^5,10^6,10^7,10^8)
                        )

fit30$mean <- fit30$num/fit30$denom
names(fit30)[3]<- "ConcMax"

plotDat +
  geom_point(data=fit30,aes(x=ConcMax,y=mean),col="red")