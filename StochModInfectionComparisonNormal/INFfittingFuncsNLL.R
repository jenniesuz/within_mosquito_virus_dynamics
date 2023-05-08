
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
                   , dat=competenceDat
                   , nSimulations=30) {
  parms <- subsParms(fit.params, fixed.params)
  nll.binom(parms=parms, dat = dat, nSims=nSimulations)                                                           ## then call likelihood function
}
#***************************************************************************************

#***********Likelihood function**********************************
nll.binom <- function(parms=virus_params()
                      ,dat=competenceDat
                      ,nSims=30){ 
  
  datAeg <- dat[dat$Moz %in% "Ae. aegypti",]
  datAlb <- dat[dat$Moz %in% "Ae. albopictus",]
  parmsAeg <- c(parms$muV 
                     ,parms$infRate1 
                     ,parms$prodRate1           
                     ,parms$cellSpread1        
                     ,parms$escapeRate1          
                     ,parms$cMax)
  parmsAlb <- c(parms$muV 
                ,parms$infRate2 
                ,parms$prodRate1              
                ,parms$cellSpread1      
                ,parms$escapeRate1      
                ,parms$cMax)
  
  #***********Aedes aegypti**************
    # replicate experiments across virus concentrations - for each virus concentration simulate 30 mosquitoes (1 experiment) 100 times
    cl <- makeCluster(detectCores()-1)
    clusterEvalQ(cl, {library(adaptivetau)})
    environment(repeatInfModel) <- .GlobalEnv
    clusterExport(cl, varlist=c("nSims","infectionModel","repeatInfModel","datAeg","parmsAeg"),
                  envir=environment())
    
    simsAeg <- parLapply(cl,1:nSims,function(y){
      #set.seed(y)
      s <- repeatInfModel(x=parmsAeg,startingVirus=10^datAeg$ConcMax)
      s$run <- y
      names(s) <-c("num","denom","conc","run")
      return(s)
    })
  stopCluster(cl)
  
    
    simsAeg <- do.call(rbind,simsAeg) 
    
    stats <- lapply(unique(simsAeg$run),function(z){
      temp <- simsAeg[simsAeg$run %in% z,]
      mod <- glm(temp$num/temp$denom~temp$conc,family="binomial",weights=temp$denom)
      return(c(temp$run[1],as.vector(coef(mod))))
    })
    
    statsAeg <- do.call(rbind.data.frame,stats)
    names(statsAeg) <- c("run","par1","par2")
    
    #*******************************
    
    #***********Aedes albopictus**************
    # replicate experiments across virus concentrations - for each virus concentration simulate 30 mosquitoes (1 experiment) 100 times
    cl <- makeCluster(detectCores()-1)
    clusterEvalQ(cl, {library(adaptivetau)})
    environment(repeatInfModel) <- .GlobalEnv
    clusterExport(cl, varlist=c("nSims","infectionModel","repeatInfModel","datAlb","parmsAlb"),
                  envir=environment())
    
    simsAlb <- parLapply(cl,1:nSims,function(y){
      #set.seed(y)
      s <- repeatInfModel(x=parmsAlb,startingVirus=10^datAlb$ConcMax)
      s$run <- y
      names(s) <-c("num","denom","conc","run")
      return(s)
    })
    stopCluster(cl)
    
    simsAlb <- do.call(rbind,simsAlb) 
    
    stats <- lapply(unique(simsAlb$run),function(z){
      temp <- simsAlb[simsAlb$run %in% z,]
      mod <- glm(temp$num/temp$denom~temp$conc,family="binomial",weights=temp$denom)
      return(c(temp$run[1],as.vector(coef(mod))))
    })
    
    statsAlb <- do.call(rbind.data.frame,stats)
    names(statsAlb) <- c("run","par1","par2")

    #***********************data************************
    modDatAeg <- glm(datAeg$NumInf/datAeg$ITotal~datAeg$ConcMax,family="binomial",weights=datAeg$ITotal)
    modDatAlb <- glm(datAlb$NumInf/datAlb$ITotal~datAlb$ConcMax,family="binomial",weights=datAlb$ITotal)

    #***************************************************
    
    llAeg <- dmvnorm(c(coef(modDatAeg))
                  ,mean=c(mean(statsAeg[,2])
                          ,mean(statsAeg[,3]))
                  ,sigma=cov(statsAeg[,2:3]),log=T)

    
    llAlb <- dmvnorm(c(coef(modDatAlb))
                     ,mean=c(mean(statsAlb[,2])
                             ,mean(statsAlb[,3]))
                     ,sigma=cov(statsAlb[,2:3]),log=T)
    
     ll <- sum(llAeg,llAlb)
  return(-ll)
}
#**********************************************************



