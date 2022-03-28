require(bbmle)  # required packages
library(parallel)
require(sensitivity)
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
  
  if (parms$muV > 1 | parms$infRate > 1  | parms$prodRate > 1000  | parms$cellSpread > 1 | parms$escapeRate > 1 ) {                                
    ll <- -1000000000
  }else{                                                      
    
 
     doseSim <- mclapply(10^dat$Conc.Min, sim.func,parms=parms)
      doseSim2 <- do.call(rbind.data.frame,doseSim)
    
    infDat <- ddply(doseSim2,.(conc,run),summarise,sumRuns=sum(inf))
    infDat$inf <- 0
    infDat$inf[infDat$sumRuns>0] <- 1
    infDat <- ddply(infDat,.(conc),summarise,prop=sum(inf)/max(run))
   
    samples <- ddply(dat,.(Conc.Min),summarise, NumInf=sum(NumInf),Denom=sum(ITotal))
    
    infDat$prop[infDat$prop %in% 0]<- 0.000000000001
    infDat$prop[infDat$prop %in% 1]<- 0.999999999999
    dbinoms <- dbinom(samples$Denom-samples$NumInf,samples$Denom,prob=1-infDat$prop,log=T)
    
    ll <- sum(dbinoms)  # log likelihood assuming data are Poisson distributed
  }

return(-ll)
}
#**********************************************************
# summarise fitted model output
modOutFunc <- function(modelOutput=doseSimAegHND){
infDat <- ddply(modelOutput,.(conc,run),summarise,sumRuns=sum(inf))
infDat$inf <- 0
infDat$inf[infDat$sumRuns>0] <- 1
infDat <- ddply(infDat,.(conc),summarise,prop=sum(inf)/max(run))
names(infDat) <- c("Conc.Min","meanInf")
return(infDat)
}
