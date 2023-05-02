
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
                   , dat=testDat
                   , nSimulations=30) {
  parms <- subsParms(fit.params, fixed.params)
  nll.binom(parms, dat = dat, nSims=nSimulations)                                                           ## then call likelihood function
}
#***************************************************************************************

#***********Likelihood function**********************************
nll.binom <- function(parms=virus_params()
                      ,dat=testDat
                      ,nSims=30){ 
  
  if (parms$muV > 1 
      | parms$infRate > 1  
      | parms$prodRate > 1000   
      | parms$cellSpread > 1
      | parms$escapeRate > 1  ) {                                
    ll <- -1000000000
  }else{                                                      
 
    # replicate experiments across virus concentrations - for each virus concentration simulate 30 mosquitoes (1 experiment) 100 times
    cl <- makeCluster(detectCores()-1)
    clusterEvalQ(cl, {library(adaptivetau)})
    environment(repeatInfModel) <- .GlobalEnv
    clusterExport(cl, varlist=c("nSims","infectionModel","repeatInfModel","dat","parms"),
                  envir=environment())
    
    sims <- parLapply(cl,1:nSims,function(y){
      #set.seed(y)
      s <- repeatInfModel(x=parms,startingVirus=10^dat$ConcMax)
      return(s)
    })
  stopCluster(cl)
    
    sims <- do.call(rbind,sims) 
    
    stats <- lapply(unique(sims$conc),function(z){
      temp <- sims[sims$conc %in% z,]
      indPs <- temp$num/temp$denom
      meanP <- sum(temp$num)/sum(temp$denom)
      betaEst <- tryCatch( { ebeta(indPs) 
      }
      ,
      error=function(cond) {
        return(c(NA,NA))
      }
      )
      if(is.na(betaEst[1])==F){
        betaEst <- as.numeric(betaEst$parameters)
      }    
      return(c(temp$conc[1],betaEst,meanP))
    })
    
    stats <- do.call(rbind.data.frame,stats)
    names(stats) <- c("ConcMax","shape1","shape2","mean")
    
    #*******************************
    
    #***********************data************************
     samples <- merge(stats,dat,by.x=c("ConcMax"))
    #***************************************************
    betaBinDat <- samples[!is.na(samples$shape1),]
    binomDat <- samples[is.na(samples$shape1),]
    
    if(length(binomDat$ConcMax>0)){
      #dBinoms <- dbinom(binomDat$NumInf,binomDat$ITotal,prob=binomDat$mean,log=T)
      ll <- -1000000000
    }else{
    
    dBetaBinoms <- dbetabinom(betaBinDat$NumInf
                              ,shape1=betaBinDat$shape1
                              ,shape2=betaBinDat$shape2
                              ,size=nSims 
                              ,log=T)
     
    
      ll <- sum(dBetaBinoms)
    }
  }
    
  return(-ll)
}
#**********************************************************



