
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
      s <- repeatInfModel(x=parmsAeg,startingVirus=10^dat$ConcMax[dat$Moz %in% "Ae. aegypti"])
      return(s)
    })
  stopCluster(cl)
    
    aeg <- do.call(rbind,aeg) 
    # establish how many runs failed
    #fails <- length(aeg$num[aeg$num %in% NA])
    aeg <- aeg[!aeg$conc %in% NA, ]
    stats <- lapply(unique(aeg$conc),function(z){
      temp <- aeg[aeg$conc %in% z,]
      indPs <- 1 - temp$num/temp$denom
      meanP <- 1 - sum(temp$num)/sum(temp$denom)
      betaEst <- tryCatch( { ebeta(indPs) 
      }
      ,
      error=function(cond) {
        message(cond)
        return(c(NA,NA))
      }
      )
      if(is.na(betaEst[1])==F){
        betaEst <- as.numeric(betaEst$parameters)
      }    
      return(c(temp$conc[1],betaEst,meanP))
    })
    
    stats <- do.call(rbind.data.frame,stats)
    names(stats) <- c("Conc.Min","shape1","shape2","mean")
    
    if(length(stats$shape1[is.na(stats$shape1)])>0){
      ll <- -1000000000
    }else{
    
    #***********************data************************
    samples <- ddply(dat,.(Moz,Conc.Min),summarise, NumInf=sum(NumInf),Denom=sum(ITotal))
    samplesAeg <- samples[samples$Moz %in% "Ae. aegypti",]
    samplesAeg <- merge(samplesAeg,stats,by.x="Conc.Min")
    #***************************************************
  
    dBetaBinomsAeg <- dbetabinom(samplesAeg$Denom-samplesAeg$NumInf
                              ,shape1=stats$shape1
                              ,shape2=stats$shape2
                              ,size=100 
                              ,log=T)
     
    ll <- sum(dBetaBinoms,na.rm=T)
    }
    }
  
  return(-ll)
}
#**********************************************************



