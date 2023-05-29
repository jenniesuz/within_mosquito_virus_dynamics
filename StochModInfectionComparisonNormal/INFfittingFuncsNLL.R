
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
      temp <- temp[!temp$num %in% NA,]
      if(length(temp$num)>0){
      mod <- glm(temp$num/temp$denom~temp$conc,family="binomial",weights=temp$denom)
      coefs <- as.vector(coef(mod))
      return(cbind.data.frame(run=temp$run[1],par1=coefs[1],par2=coefs[2]))
      }
      else{return(cbind.data.frame(run=temp$run[1],par1=NA,par2=NA))}
    })
    
    statsAeg <- do.call(rbind.data.frame,stats)
    
    
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
      temp <- temp[!temp$num %in% NA,]
      if(length(temp$num)>0){
        mod <- glm(temp$num/temp$denom~temp$conc,family="binomial",weights=temp$denom)
        coefs <- as.vector(coef(mod))
        return(cbind.data.frame(run=temp$run[1],par1=coefs[1],par2=coefs[2]))
      }
      else{return(cbind.data.frame(run=temp$run[1],par1=NA,par2=NA))}
    })
    
    statsAlb <- do.call(rbind.data.frame,stats)

    #***********************data************************
    modDatAeg <- glm(datAeg$NumInf/datAeg$ITotal~datAeg$ConcMax,family="binomial",weights=datAeg$ITotal)
    modDatAlb <- glm(datAlb$NumInf/datAlb$ITotal~datAlb$ConcMax,family="binomial",weights=datAlb$ITotal)

    #***************************************************
    statsAeg<- statsAeg[!statsAeg$par1 %in% NA,]
    statsAlb<- statsAlb[!statsAlb$par1 %in% NA,]
    
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






#**********************for mle2******************
#***********Likelihood function**********************************
nll.binom.mle2 <- function(logmuV=log(0.1)
                      ,loginfRate=log(10^-8)
                      ,prodRate=1
                      ,cellSpread=10^-4
                      ,escapeRate=0.05
                      ,cMax=400
                      ,dat=competenceDat
                      ,nSims=30){ 
  if(exp(logmuV)>1|exp(loginfRate)>0.005){ll<- -10000000}else{
  parms <- c("muV"=exp(logmuV) 
                ,"infRate"=exp(loginfRate) 
                ,"prodRate"=prodRate          
                ,"cellSpread"=cellSpread       
                ,"escapeRate"=escapeRate         
                ,"cMax"=cMax)

  # replicate experiments across virus concentrations - for each virus concentration simulate 30 mosquitoes (1 experiment) 100 times
  cl <- makeCluster(detectCores()-1)
  clusterEvalQ(cl, {library(adaptivetau)})
  environment(repeatInfModel) <- .GlobalEnv
  clusterExport(cl, varlist=c("nSims","infectionModel","repeatInfModel","dat","parms"),
                envir=environment())
  
  sims <- parLapply(cl,1:nSims,function(y){
    #set.seed(y)
    s <- repeatInfModel(x=parms,startingVirus=10^dat$ConcMax)
    s$run <- y
    names(s) <-c("num","denom","conc","run")
    return(s)
  })
  stopCluster(cl)
  
  sims <- do.call(rbind,sims) 
  
  stats <- lapply(unique(sims$run),function(z){
    temp <- sims[sims$run %in% z,]
    temp <- temp[!temp$num %in% NA,]
    if(length(temp$num)>0){
      mod <- glm(temp$num/temp$denom~temp$conc,family="binomial",weights=temp$denom)
      coefs <- as.vector(coef(mod))
      return(cbind.data.frame(run=temp$run[1],par1=coefs[1],par2=coefs[2]))
    }
    else{return(cbind.data.frame(run=temp$run[1],par1=NA,par2=NA))}
  })
  
  stats <- do.call(rbind.data.frame,stats)
  
  #***********************data************************
  modDat <- glm(dat$NumInf/dat$ITotal~dat$ConcMax,family="binomial",weights=dat$ITotal)

  #***************************************************
  stats <- stats[!stats$par1 %in% NA,]
  if(length(stats)==0){ll<- -10000000000}else{
  ll <- dmvnorm(c(coef(modDat))
                   ,mean=c(mean(stats[,2])
                           ,mean(stats[,3]))
                   ,sigma=cov(stats[,2:3]),log=T)
  }
}

  return(-ll)

}
#**********************************************************

