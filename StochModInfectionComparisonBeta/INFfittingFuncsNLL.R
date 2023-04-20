
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



#************************parameters for two 'treatments'****************************
virus_params <- function(   muV = 0.1
                            ,infRate1 = 10^-7.5
                            ,infRate2 = 10^-7.5
                            ,prodRate1 =   1    
                            ,prodRate2 = 1
                            ,cellSpread1 = 10^-4  
                            ,cellSpread2 = 10^-4
                            ,escapeRate1 = 0.05
                            ,escapeRate2 = 0.05
                            ,cMax = 400  
)
  return(as.list(environment()))
#*****************************************************************************

#**************************Function to wrap fitting and simulations for iteration*************
#repeatModelRunFunc <- function(x
#                               ,conc=dat$Conc.Min
#                               ,p1=parms$muV
#                               ,p2=parms$infRate1
#                               ,p3=parms$prodRate1
#                               ,p4=parms$cellSpread1
#                               ,p5=parms$escapeRate1
#                               ,p6=parms$cMax
#                               ){
#  a<-x
#  doseSim <- lapply(10^conc, sim.func
#                    ,muV=p1
#                    ,infRate=p2
#                    ,prodRate=p3
#                    ,cellSpread=p4
#                    ,escapeRate=p5
#                    ,cMax=p6)
#  doseSim2 <- do.call(rbind.data.frame,doseSim)
#  infDat <- ddply(doseSim2,.(conc,run),summarise,sumRuns=sum(inf))
#  infDat$inf <- 0
#  infDat$inf[infDat$sumRuns>0] <- 1
#  infDat$run2 <- x
#  return(infDat)
#}
#********************************************************************************


#***********Likelihood function**********************************
nll.binom <- function(parms=virus_params()
                      ,dat=competenceDat){ 
  
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
 
    # separate parameters for each mosquito
    parmsAeg  <- c("muV"=parms$muV,"infRate"=parms$infRate1,"prodRate"=parms$prodRate1,"cellSpread"=parms$cellSpread1,"escapeRate"=parms$escapeRate1,"cMax"=parms$cMax)
    parmsAlb <- c("muV"=parms$muV,"infRate"=parms$infRate2,"prodRate"=parms$prodRate2,"cellSpread"=parms$cellSpread2,"escapeRate"=parms$escapeRate2,"cMax"=parms$cMax)
    
    #********aegypti******
    # replicate experiments across virus concentrations - for each virus concentration simulate 30 mosquitoes (1 experiment) 100 times
    cl <- makeCluster(detectCores()-1)
    clusterEvalQ(cl, {library(adaptivetau)})
    environment(repeatInfModel) <- .GlobalEnv
    clusterExport(cl, varlist=c("infectionModel","repeatInfModel","dat","parmsAeg"),
                  envir=environment())
    
    aeg <- parLapply(cl,1:100,function(y){
      #set.seed(y)
      s <- repeatInfModel(x=parms,startingVirus=10^dat$ConcMax[dat$Moz %in% "Ae. aegypti"])
      return(s)
    })
  stopCluster(cl)
    
    aeg <- do.call(rbind,aeg) 
    # establish how many runs failed
    #fails <- length(aeg$num[aeg$num %in% NA])

    stats <- lapply(unique(aeg$conc),function(z){
      temp <- aeg[aeg$conc %in% z,]
      numM <- mean(temp$num,na.rm=T)
    })
    
    #*data*
    samples <- ddply(dat,.(Moz,Conc.Min),summarise, NumInf=sum(NumInf),Denom=sum(ITotal))
    samplesAeg <- samples[samples$Moz %in% "Ae. aegypti",]
    # samplesAlb <- samples[samples$Moz %in% "Ae. albopictus",]
    samples <- rbind.data.frame(samplesAeg,samplesAlb)
    
    infDat[infDat %in% 0]<- 0.000000000001
    infDat[infDat%in% 1]<- 0.999999999999
    dbinoms <- dbinom(samples$Denom-samples$NumInf,samples$Denom,prob=1-infDat,log=T)
    
    ll <- sum(dbinoms)  # log likelihood assuming data are Poisson distributed
    }
  }
return(-ll)
}
#**********************************************************



