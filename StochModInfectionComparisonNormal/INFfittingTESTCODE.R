library(here)
library(binom)
library(bbmle)
library(ggplot2)
library(adaptivetau)
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


#***************explore likelihoods************
objFXN(fit.params=init.pars                                                          ## paramters to fit
       , fixed.params =virus_params()                                      ## fixed paramters
       , dat=testDat
       , nSimulations=30)


# calculate likelihood for range of values 
# have to constrain the data to be between 1 infected and n-1


start <- Sys.time()
testLik30 <- lapply(10^-seq(7,8,0.05),function(x){
  loglik <- nll.binom(parms=virus_params(muV=0.1
                               ,infRate=x)
            ,dat=testDat
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
# There is a trade-off between number of simulations used to
# estiamte the beta parameters and therefore accuracy of those
# estimations and the time it takes to fit the model


#********initial parameter values****
init.pars.fit <- c(
  log_infRate=log(10^-7) 
 # ,log_muV=log(0.015)
)


# tol 0.0001 too small
# tol 0.001 took 10 mins to 49 mins

#**************************************
start <- Sys.time()

#********optimise*******
optim.vals <- optim(par = init.pars.fit
                    , objFXN
                    , fixed.params = virus_params()
                    , nSimulations = 30
                    , dat = testDat
                    , control = list(trace = 3
                                     ,abstol=0.001
                                     ,reltol=0.001
                                     )
                    , method = "Nelder-Mead" # 
                    , hessian = F)
  
 end <- Sys.time()
end-start






# library(lme4)
# start <- Sys.time()
# 
# #********optimise*******
# optim.vals2 <- NelderMead(, objFXN
#                      , x0 = init.pars.fit
#                      , fixed.params = virus_params()
#                      , lower=c(log(10^-8))
#                      , upper=c(log(10^-6))
#                      , xst=0.02
#                      , nSimulations = 30
#                      , dat = testDat
#                      , control = list(trace = 3
#                                       ,FtolAbs=0.001
#                                       ,FtolRel=0.001
#                      )
#                      )
# 
# end <- Sys.time()
# end-start



saveRDS(optim.vals,"OptimVals30sims100mozInfRateOnly.RDS") # when 100 maxit takes c. 20-30 mins # when 200 sims took 1 hour

fisherInfMatrix <- solve(optim.vals$hessian) ## invert the Hessian, to estimate the covar-var matrix of parameter estimates
fisherInfMatrix

# Finds the critical z value
conf.level <- 0.95
crit <- qnorm((1 + conf.level)/2)

ci <- optim.vals$par[1] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[1, 1]))
log10(exp(ci))


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


saveRDS(optim.vals,"OptimVals30sims30mozGenSA.RDS") # when 100 maxit takes c. 20-30 mins # when 200 sims took 1 hour

# when maxtime set to 60 sections 30 sims and 30 moz tooc 1.4

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