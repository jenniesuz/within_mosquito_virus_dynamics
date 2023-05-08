

#***********Likelihood function**********************************
nll.binom <- function(
                      ,parms=virus_params()
                      ,dat=testDat
                      ,nSims=30){ 
  
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
   # binomDat <- samples[is.na(samples$shape1),]
    
   # if(length(binomDat$ConcMax>0)){
      #dBinoms <- dbinom(binomDat$NumInf,binomDat$ITotal,prob=binomDat$mean,log=T)
  #    ll <- -1000000000
  #  }else{
    
    dBetaBinoms <- dbetabinom(betaBinDat$NumInf
                              ,shape1=betaBinDat$shape1
                              ,shape2=betaBinDat$shape2
                              ,size=nSims 
                              ,log=T)
     
    
      ll <- sum(dBetaBinoms)
   # }
  #}
    
  return(-ll)
}
#**********************************************************



