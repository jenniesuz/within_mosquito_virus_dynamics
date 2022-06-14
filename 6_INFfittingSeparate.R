library(gridExtra)
library(ggplot2)
library(binom)
library(beepr)
library(parallel)
library(bbmle)

library(GillespieSSA)

sim.func <- function(x=10000,muV=0.1,infRate=3.162278e-08,prodRate=1,cellSpread=1e-04
,escapeRate=0.05,cMax=400){
  initial <- c(Gv=x*0.003,Mci=0,Mv=0) # initial state # assume input virus concentration then mosquito imbibes c. 3ul
  
  finalTime <- 124
  
  nsims <- 30
  
  a <- c("muV*Gv"
         ,"Gv*infRate*(cMax-Mci)"
         ,"cellSpread*Mci*(cMax-Mci)"
         ,"prodRate*Mci"
         ,"muV*Mv"
         ,"escapeRate*Mv"
  ) # character vector of propensity functions
  #
  nu <- cbind(c(-1,0,0)
              ,c(-1,+1,0)
              ,c(0,+1,0)
              ,c(0,0,+1)
              ,c(0,0,-1)
              ,c(0,0,-1)
  )
  
  repSim <- lapply(1:nsims,function(z){
    out <- ssa(initial,a,nu,parms,tf=finalTime,method="ETL")
    
    dat<-data.frame(out$dat)
    dat$run <- z
    
    dat$inf <- 0
    
    if(dat$Mv[length(dat$Mv)]>0){dat$inf <- 1}
    return(dat)
  })
  
  repSims <- do.call(rbind.data.frame,repSim)
  repSims$conc <- x
  repSims$MciProp <- repSims$Mci/as.numeric(parms[6])
  return(repSims)
}



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


#***********Likelihood function**********************************
nll.binom <- function(parms=virus_params(),dat=competenceDat){ 
  
  if (parms$muV > 1 | parms$infRate > 1  | parms$infRate1 > 1  | parms$prodRate > 1000  | parms$prodRate1 > 1000  | parms$cellSpread > 1 | parms$cellSpread1 > 1 |  parms$escapeRate > 1 | parms$escapeRate1 > 1  ) {                                
    ll <- -1000000000
  }else{                                                      
    
    muV <- as.numeric(parms[1])
    cMax <- as.numeric(parms[10])
    # aeg
    infRate <- as.numeric(parms[2])
    prodRate <- as.numeric(parms[4])
    cellSpread <- as.numeric(parms[6])
    escapeRate <- as.numeric(parms[8])
   
    doseSim <- mclapply(round(10^dat$Conc.Min,0), sim.func,muV,infRate,prodRate,cellSpread,escapeRate,cMax)
    doseSim2 <- do.call(rbind.data.frame,doseSim)
    infDat <- ddply(doseSim2,.(conc,run),summarise,sumRuns=sum(inf))
    infDat$inf <- 0
    infDat$inf[infDat$sumRuns>0] <- 1
    infDatAeg <- ddply(infDat,.(conc),summarise,prop=sum(inf)/max(run))
    # alb
    doseSim <- mclapply(10^dat$Conc.Min, sim.func,parms=parms[c(1,3,5,7,9,10)])
    doseSim2 <- do.call(rbind.data.frame,doseSim)
    infDat <- ddply(doseSim2,.(conc,run),summarise,sumRuns=sum(inf))
    infDat$inf <- 0
    infDat$inf[infDat$sumRuns>0] <- 1
    infDatAlb <- ddply(infDat,.(conc),summarise,prop=sum(inf)/max(run))
    
    infDat <- rbind.data.frame(infDatAeg,infDatAlb)
    
    samples <- ddply(dat,.(Moz,Conc.Min),summarise, NumInf=sum(NumInf),Denom=sum(ITotal))
    
    
    infDatAeg$prop[infDatAeg$prop %in% 0]<- 0.000000000001
    infDatAeg$prop[infDatAeg$prop %in% 1]<- 0.999999999999
    infDatAlb$prop[infDatAlb$prop %in% 0]<- 0.000000000001
    infDatAlb$prop[infDatAlb$prop %in% 1]<- 0.999999999999
    
    dbinoms <- dbinom(samples$Denom-samples$NumInf,samples$Denom,prob=1-infDat$prop,log=T)
    
    ll <- sum(dbinoms)  # log likelihood assuming data are Poisson distributed
  }
  
  return(-ll)
}
#**********************************************************




#************************parameters****************************
virus_params <- function(   muV = 0.1
                            ,infRate = 10^-7.5
                            ,infRate1 = 10^-7.5
                            ,prodRate =   1    
                            ,prodRate1 = 1
                            ,cellSpread = 10^-4  
                            ,cellSpread1 = 10^-4
                            ,escapeRate = 0.05
                            ,escapeRate1 = 0.05
                            ,cMax = 400  
)
  return(as.list(environment()))
#*************************************************************
#******************data to fit to********************
competenceDat <- read.csv("datCiotaOnyango.csv")
competenceDat <- competenceDat[competenceDat$Ref %in% "Ciota 2017",]
bin <- binom.confint(x=competenceDat$NumInf,n=competenceDat$ITotal,method="exact")
competenceDat$meanInf <- bin$mean
competenceDat$lowerInf <- bin$lower
competenceDat$upperInf <- bin$upper
#****************************************************




fitmod <- mle2(function(par1,par2){nll.pois(par1,par2)},
               start=list(par1=log(0.124),par2=log(0.018))
               ,data=list(tsetse)) # run with with method = "SANN" first
