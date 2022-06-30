require(bbmle)  # required packages
library(parallel)
library(plyr)
#********enable multiple parameters to be fit and substituted back in parameter list****
subsParms <- function(fit.params, fixed.params=virus_params())
  
  within(fixed.params, {
    loggedParms <- names(fit.params)[grepl('log_', names(fit.params))]
    
    unloggedParms <- names(fit.params)[!grepl('log_', names(fit.params))]  
    
    for(nm in unloggedParms) assign(nm, as.numeric(fit.params[nm]))
    
    for(nm in loggedParms) assign(gsub('log_','',nm), exp(as.numeric(fit.params[nm])))
    
    rm(nm, loggedParms, unloggedParms)
  })           

## Make likelihood a function of fixed and fitted parameters.
objFXN <- function(fit.params                                                          ## paramters to fit
                   , fixed.params =virus_params()                                     ## fixed paramters
                   , dat=temps.count) {
  parms <- subsParms(fit.params, fixed.params)
  nll.binom(parms, dat = dat)                                                           ## then call likelihood function
}
#***************************************************************************************


#***********Likelihood function**********************************
nll.binom <- function(parms=virus_params(),dat=competenceDat){ 
  
  if (
        parms$prodRate > 1000  
      | parms$cellSpread > 1 
      | parms$escapeRate > 1 ) {                                
    ll <- -1000000000
  }else{     
    
      sim <- sim.func(10^dat$Conc.Min[1],parms=parms)
      
      #**********first remove simulations where an infection didn't occur**********************
      inf <- ddply(sim,.(run),summarise,occurred=sum(inf))  # establish if infection occurred
      inf$occurred[inf$occurred>0] <- 1
      noInf <- inf$run[inf$occurred == 0]                         # note run and concs where infection didn't occur
      # remove these from the dissemination dataset
      sim <- sim[!sim$run %in% noInf,]
      #*******************************************************************************
      
      tempDiss <- dissSummaryFunc(sim)
      subSim <- tempDiss[tempDiss$time %in% as.numeric(dat$DPIDissorTrans),]
      
    subSim$proportionDisseminated[subSim$proportionDisseminated %in% 0]<- 0.000000000001
    subSim$proportionDisseminated[subSim$proportionDisseminated %in% 1]<- 0.999999999999
    dbinoms <- dbinom(dat$NumInf-dat$NumDiss,dat$NumInf,prob=1-subSim$proportionDisseminated,log=T)
    
    ll <- sum(dbinoms)  # log likelihood assuming data are Poisson distributed
  }

return(-ll)
}
#**********************************************************


