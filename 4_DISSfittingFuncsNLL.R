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
    
      sim <- sim.func(10^round(dat$Conc.Min[1],0),parms=parms)
      
      sim$days <- round(sim$t/24,4)
        tempDiss <- lapply(unique(sim$days),function(z){ # for all runs sample a time point
          subDat <- sim[sim$days %in% z,]               # at that time point total number of runs that
          HciInf <- length(subDat$Hci[subDat$Hci>0])   # have Hci>0
          MvInf <- length(subDat$Mv[subDat$Mv>0])
          if(MvInf==0){propDiss<-0}else{
            propDiss <- HciInf/ MvInf 
          }#proportion 
          return(c(z,propDiss))
        })
      tempDiss <- do.call(rbind.data.frame,tempDiss)
      names(tempDiss) <- c("time","proportionDisseminated")
  
      subSim <- tempDiss[tempDiss$days %in% as.numeric(dat$DPIDissorTrans),]
      
    subSim$proportionDisseminated[subSim$proportionDisseminated %in% 0]<- 0.000000000001
    subSim$proportionDisseminated[subSim$proportionDisseminated %in% 1]<- 0.999999999999
    dbinoms <- dbinom(dat$NumInf-dat$NumDiss,dat$NumInf,prob=1-subSim$proportionDisseminated,log=T)
    
    ll <- sum(dbinoms)  # log likelihood assuming data are Poisson distributed
  }

return(-ll)
}
#**********************************************************


