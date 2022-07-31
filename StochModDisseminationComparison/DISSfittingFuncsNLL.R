
library(GillespieSSA)

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
                   , fixed.params =virus_params()                                      ## fixed paramters
                   , dat=temps.count) {
  parms <- subsParms(fit.params, fixed.params)
  nll.binom(parms, dat = dat)                                                           ## then call likelihood function
}
#***************************************************************************************


#log_muV 1.544027e-01 log_infRate1  2.338052e-09 log_infRate2  5.494484e-08 

#************************parameters for two 'treatments'****************************
virus_params <- function(   muV = 1.544027e-01
                            ,infRate1 = 2.338052e-09
                            ,infRate2 = 5.494484e-08 
                            ,prodRate1 =  50   
                            ,prodRate2 = 10
                            ,cellSpread1 = 10^-3.5 
                            ,cellSpread2 = 10^-6
                            ,escapeRate1 = 0.09
                            ,escapeRate2 = 0.005
                            ,cMax = 400  
                            ,hMax= 900
)
  return(as.list(environment()))
#*****************************************************************************

#**************************Function to wrap fitting and simulations for iteration*************
repeatModelRunFunc <- function(x
                               ,dpid=dat$DPIDissorTrans
                               ,conc=dat$Conc.Min
                               ,p1=parms$muV
                               ,p2=parms$infRate1
                               ,p3=parms$prodRate1
                               ,p4=parms$cellSpread1
                               ,p5=parms$escapeRate1
                               ,p6=parms$cMax
                               ,p7=parms$hMax
                               ){
  a<-x
  sim <- sim.func(x=10^conc[1]
                    ,muV=p1
                    ,infRate=p2
                    ,prodRate=p3
                    ,cellSpread=p4
                    ,escapeRate=p5
                    ,cMax=p6
                    ,hMax=p7)

  #**********first remove simulations where an infection didn't occur**********************
  inf <- ddply(sim,.(run),summarise,occurred=sum(inf))  # establish if infection occurred
  inf$occurred[inf$occurred>0] <- 1
  noInf <- inf$run[inf$occurred == 0]                         # note run and concs where infection didn't occur
  # remove these from the dissemination dataset
  sim <- sim[!sim$run %in% noInf,]
  #*******************************************************************************
  tempDiss <- dissSummaryFunc(sim)

  subSim <- tempDiss[tempDiss$time %in% as.numeric(dpid),]
  subSim$run2 <- x
  return(subSim)
}
#********************************************************************************


#***********Likelihood function**********************************
nll.binom <- function(parms=virus_params(),dat=competenceDat){ 
  
  if (parms$muV > 1 
      | parms$infRate1 > 1  
      | parms$infRate2 > 1  
      | parms$prodRate1 > 1000  
      | parms$prodRate2 > 1000  
      | parms$cellSpread1 > 1
      | parms$cellSpread2 > 1
      | parms$escapeRate1 > 1 
      | parms$escapeRate2 > 1 ) {                                
    ll <- -1000000000
  }else{                                                      
 
     # aeg
     sim <- mclapply(1:30,repeatModelRunFunc    # 30 simulations to find average given set of parameter values
                           ,p2=parms$infRate1
                           ,dpid=dat$DPIDissorTrans[dat$Moz %in% "Ae. aegypti"]
                           ,conc=dat$Conc.Min[dat$Moz %in% "Ae. aegypti"]
                        ,p1=parms$muV
                        ,p3=parms$prodRate1
                        ,p4=parms$cellSpread1
                        ,p5=parms$escapeRate1
                        ,p6=parms$cMax
                        ,p7=parms$hMax)
   
     

    sim2 <- do.call(rbind.data.frame,sim)
    sim3 <- ddply(sim2,.(run2,time),summarise,meanDiss=mean(HciInf),propDiss=mean(HciInf)/SampleSize)
    aegSim <- sim3
    
    #alb
    sim <- mclapply(1:30,repeatModelRunFunc    # 30 simulations to find average given set of parameter values
                    ,p2=parms$infRate2
                    ,dpid=dat$DPIDissorTrans[dat$Moz %in% "Ae. albopictus"]
                    ,conc=dat$Conc.Min[dat$Moz %in% "Ae. albopictus"]
                    ,p1=parms$muV
                    ,p3=parms$prodRate2
                    ,p4=parms$cellSpread2
                    ,p5=parms$escapeRate2
                    ,p6=parms$cMax
                    ,p7=parms$hMax)
    
    
    
    sim2 <- do.call(rbind.data.frame,sim)
    sim3 <- ddply(sim2,.(run2,time),summarise,meanDiss=mean(HciInf),propDiss=mean(HciInf)/SampleSize)
    albSim <- sim3
    
    aegSim$Moz <- "Ae. aegypti"
    albSim$Moz <- "Ae. albopictus"
    sims <- rbind.data.frame(aegSim,albSim)
    
    dat <- rbind.data.frame(dat[dat$Moz %in% "Ae. aegypti",],dat[dat$Moz %in% "Ae. albopictus",])
    
    sims$propDiss[sims$propDiss %in% 0]<- 0.000000000001
    sims$propDiss[sims$propDiss %in% 1]<- 0.999999999999
    dbinoms <- dbinom(dat$NumInf-dat$NumDiss,dat$NumInf,prob=1-sims$propDiss,log=T)
    
    ll <- sum(dbinoms)
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
