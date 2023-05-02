library(here)
library(binom)
library(bbmle)
library(emdbook)
library(EnvStats)
library(optimx)
library(ggplot2)
#****model code****
source(here("StochModInfectionComparisonBeta//INFmodelFunc.R"))
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
                              ,infRate=10^-7.4)
                              ,dat=testDat
                              ,nSims=30) 

# calculate likelihood for range of values - need to avoid
# values where no simulations result in infection for a given concentration
# or all simulations result in infection at a given concentration
# but where data observe infection occurring - these are out of bounds
# so these need to be set at the start of the fitting process as input parameters
# and from consultation of the data

# have to constrain the data to be between 1 infected and n-1
# have to constrain the model parameters to the space in which the 
# simulations follow a beta distribution not binomial. - if can't calculate 
# a beta distribution is it better to reutrn a very low likelihood
# as saying that the parameter values are out of bounds - if the
# model is saying something is 100% certain then can't be right
# parameterisation....therefore it becomes important to do more
# 'experiments' or more 'mosquitoes' to ensure that our
# near-zero and near-one but not quite either are represented.
# however this comes at a computational cost... can we still 
# get there if we use lower sims and just put a restriction
# on the liklihood?

start <- Sys.time()
testLik30 <- lapply(c(10^-5
                      ,10^-5.5
                      ,10^-6
                      ,10^-6.5
                      ,10^-7
                      ,10^-7.5
                      ,10^-8
                      ,10^-8.5
                      ,10^-9),function(x){
  loglik <- nll.binom(parms=virus_params(muV=0.1
                               ,infRate=x)
            ,dat=testDat
            ,nSims=30) 
  return(loglik)
})

testLik30 <- do.call(c,testLik30)
end <- Sys.time()
end-start

# 9 parameter values and 30 simulations takes 1.9mins

# There is a trade-off between number of simulations used to
# estiamte the beta parameters and therefore accuracy of those
# estimations and the time it takes to fit the model
plot( c(5
  ,5.5
  ,6
  ,6.5
  ,7
  ,7.5
  ,8
  ,8.5
  ,9), log10(testLik30))

# could do this to get an initial idea of bounds then repeat to narrow
# in

start <- Sys.time()
testLik30 <- lapply(10^seq(-9,-6,0.05),function(x){
                        loglik <- nll.binom(parms=virus_params(muV=0.1
                                                               ,infRate=x)
                                            ,dat=testDat
                                            ,nSims=30) 
                        return(loglik)
                      })

testLik30 <- do.call(c,testLik30)
end <- Sys.time()
end-start


liklihoods <- cbind.data.frame("infRate"=10^seq(-9,-6,0.05),"ll"=testLik30)
saveRDS(liklihoods,"likelihoodWhenTrueInfNeg7pt5.RDS")
#What happens when repeat above with more simulated experiments?
# goes from c. 2 mins to c. 3-4 mins so not too bad BUT 
# interestingly has produced the same results of high...why?!?!

ggplot(liklihoods) +
  geom_line(aes(log10(infRate),ll)) +
  ylim(0,100)
















#********initial parameter values****
init.pars.fit <- c(
  log_infRate=log(10^-9)  
)
  
start <- Sys.time()

#********optimise*******
optim.vals <- optimx(par = init.pars.fit
                    , objFXN
                    , fixed.params = virus_params()
                    , nSimulations = 100
                    , dat = testDat
                    , control = list(trace = 3
                                     )
                    , method = "Nelder-Mead" # 
                    , hessian = T)
  
 end <- Sys.time()
end-start

saveRDS(optim.vals,"OptimVals30sims100mozInfRateOnly.RDS") # when 100 maxit takes c. 20-30 mins # when 200 sims took 1 hour


# conv code 10 for 30 mosquitoes and 100 mosquitoes and 30 simulated experiments per dose takes 2 hours with 100 moz
# try for 30 mosquitoes adn 100 sims

log10(exp(coef(optim.vals))[1])
log10(exp(coef(optim.vals))[2])


#The Hessian matrix gives you the curvature of the likelihood function at the maximum likelihood
## estimate (MLE) of the fitted paramters. In other words, it tells you the second derivative around
## MLE, which can be used to estimate the covariance variance matrix of the MLE. This estimate of
## the covariance varance matrix is known as the Fisher information matrix and can be obtained by
## inverting the Hessian.
fisherInfMatrix <- solve(optim.vals$hessian) ## invert the Hessian, to estimate the covar-var matrix of parameter estimates


#***********try GenSA

library(GenSA)#

init.pars.fit <- c(
  log_infRate=log(10^-7)  
)

start <- Sys.time()
optim.valsGSA <- GenSA(par = init.pars.fit
                       ,fn = objFXN
                       ,dat=testDat
                       ,nSimulations = 30
                       ,fixed.params = virus_params()
                       ,lower=c(log(10^-15))
                       ,upper=c(log(10^-5))
                       ,control=list(max.time=60,smooth=F) #,thresold.stop=1)
                       
                      )


end <- Sys.time()
end-start

saveRDS(optim.vals,"OptimVals30sims30mozGenSA.RDS") # when 100 maxit takes c. 20-30 mins # when 200 sims took 1 hour

# when maxtime set to 60 sections 30 sims and 30 moz tooc 1.3 hours
# giving a strange minimum of 0 makes no sense

#********************************

library(optimization)

start <- Sys.time()
optim.valsSA <- optim_sa(fun=objFXN
                         ,start=init.pars.fit
                         ,maximization = FALSE
                         ,trace = TRUE
                         ,lower=c(10^-15)
                         ,upper=c(10^-5)
                        # ,dat=testDat,nSimulations=30,fixed.params=virus_params()
                         )


end <- Sys.time()
end-start







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