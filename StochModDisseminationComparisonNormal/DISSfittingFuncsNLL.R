
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
                   ,nSimulations=30) {
  parms <- subsParms(fit.params, fixed.params)
  nll.binom(parms=parms, dat = dat, nSims=nSimulations)                                                           ## then call likelihood function
}
#***************************************************************************************
#log_muV 1.544027e-01 log_infRate1  2.338052e-09 log_infRate2  5.494484e-08 

#************************parameters for two 'treatments'****************************
virus_params <- function(   muV = 0.1
                            ,infRate1 = 10^-8 # 2.338052e-09
                            ,infRate2 = 5.494484e-08 
                            ,prodRate1 = 25  # 50   
                            ,prodRate2 = 10
                            ,cellSpread1 = 0.00005 # 10^-3.5 
                            ,cellSpread2 = 10^-6
                            ,escapeRate1 = 0.005 # 0.0005
                            ,escapeRate2 = 0.0005
                            ,cMax = 400  
                            ,hMax= 900
)
  return(as.list(environment()))
#*****************************************************************************

#***********Likelihood function**********************************
nll.binom <- function(parms=virus_params()
                      ,dat=competenceDat
                      ,nSims=30){ 
  
  if (parms$prodRate1 > 1000  
      | parms$prodRate2 > 1000  
      | parms$cellSpread1 > 1
      | parms$cellSpread2 > 1
      | parms$escapeRate1 > 1 
      | parms$escapeRate2 > 1 ) {                                
    ll <- -1000000000
  }else{                                                      
 
    datAeg <- dat[dat$Moz %in% "Ae. aegypti",]
    datAlb <- dat[dat$Moz %in% "Ae. albopictus",]
    parmsAeg <- c(parms$muV 
                  ,parms$infRate1 
                  ,parms$prodRate1           
                  ,parms$cellSpread1        
                  ,parms$escapeRate1          
                  ,parms$cMax
                  ,parms$hMax)
    parmsAlb <- c(parms$muV 
                  ,parms$infRate1
                  ,parms$prodRate1              
                  ,parms$cellSpread1      
                  ,parms$escapeRate1      
                  ,parms$cMax
                  ,parms$hMax)
    
    
    #***********Aedes aegypti**************
    # replicate experiments across virus concentrations - for each virus concentration simulate 30 mosquitoes (1 experiment) 100 times
    cl <- makeCluster(detectCores()-1)
    clusterEvalQ(cl, {library(adaptivetau)})
    environment(repeatModel) <- .GlobalEnv
    clusterExport(cl, varlist=c("nSims","infDissModel","dissSummaryFunc" ,"repeatModel","datAeg","parmsAeg"),
                  envir=environment())
    
    simsAeg <- parLapply(cl,1:nSims,function(y){
      #set.seed(y)
      s <- repeatModel(x=parmsAeg,startingVirus=10^datAeg$ConcMax[1])
      s <- do.call(rbind,s)
      # get rid of runs where midgut infection didn't occur
      s <- s[s$inf %in% 1,]   # make sure calculating the probability of dissemination given infection
      calcDisseminations <- dissSummaryFunc(s)
      calcDisseminations$run <- y
      return(calcDisseminations)
    })
    stopCluster(cl)
    
    simsAeg <- do.call(rbind,simsAeg) 
   
     print(length(unique(simsAeg$run)))
    
    stats <- lapply(unique(simsAeg$run),function(z){
      temp <- simsAeg[simsAeg$run %in% z,]
      temp <- temp[!temp$num %in% NA,]
      # cut of when all disseminated
      temp <- temp[temp$numberRunsDisseminated<temp$totalSize,]
      if(length(temp$num)>0){
        mod <- glm(temp$numberRunsDisseminated/temp$totalSize~temp$time,family="binomial",weights=temp$totalSize)
        coefs <- as.vector(coef(mod))
        return(cbind.data.frame(run=temp$run[1],par1=coefs[1],par2=coefs[2]))
      }
      else{return(cbind.data.frame(run=temp$run[1],par1=NA,par2=NA))}
    })
    
    statsAeg <- do.call(rbind.data.frame,stats)
    
    print(length(statsAeg[,1]))
    
    #********observed data****************
    modDatAeg <- glm(datAeg$NumDiss/datAeg$NumInf~datAeg$DPIDissorTrans,family="binomial",weights=datAeg$NumInf)
    modDatAlb <- glm(datAlb$NumDiss/datAlb$NumInf~datAlb$DPIDissorTrans,family="binomial",weights=datAlb$NumInf)
    
    statsAeg <- statsAeg[!statsAeg$par1 %in% NA,]
    statsAlb <- statsAlb[!statsAlb$par1 %in% NA,]
    if(length(statsAeg)==0|length(statsAlb==0)){ll<- -10000000000}else{
      llAeg <- dmvnorm(c(coef(modDatAeg))
                       ,mean=c(mean(statsAeg[,2])
                               ,mean(statsAeg[,3]))
                       ,sigma=cov(statsAeg[,2:3]),log=T)
      
      
      llAlb <- dmvnorm(c(coef(modDatAlb))
                       ,mean=c(mean(statsAlb[,2])
                               ,mean(statsAlb[,3]))
                       ,sigma=cov(statsAlb[,2:3]),log=T)
      
      ll <- sum(llAeg,llAlb)
    }

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
