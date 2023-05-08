
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
                      ,dat=competenceDat
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
      s$run <- y
      return(s)
    })
  stopCluster(cl)
    
    sims <- do.call(rbind,sims) 
    
    stats <- lapply(unique(sims$run),function(z){
      temp <- sims[sims$run %in% z,]
      mod <- glm(temp$num/temp$denom~temp$conc,family="binomial",weights=temp$denom)
      return(c(temp$run[1],as.vector(coef(mod))))
    })
    
    stats <- do.call(rbind.data.frame,stats)
    names(stats) <- c("run","par1","par2")
    
    #*******************************
    
    #***********************data************************
     modDat <- glm(dat$NumInf/dat$ITotal~dat$ConcMax,family="binomial",weights=dat$ITotal)
    #***************************************************
 
    ll <- dmvnorm(c(coef(modDat))
            ,mean=c(mean(stats[,2])
                    ,mean(stats[,3]))
            ,sigma=cov(stats[,2:3]),log=T)
    
    
  return(-ll)
}
#**********************************************************



