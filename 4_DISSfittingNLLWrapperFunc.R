library(gridExtra)
library(ggplot2)
library(binom)
library(beepr)

source("4_DISSfittingFuncsNLL.R")
source("4_DISSmodelFunc.R")

virus_params <- function(   muV=0.1
                            ,infRate=10^-7   
                            ,prodRate =   10 
                            ,cellSpread = 10^-3.5  
                            ,escapeRate = 0.05   
                            ,cMax = 400            
                            ,hMax = 900             
)
  return(as.list(environment()))


#******************data to fit to********************
competenceDat <- read.csv("datCiotaOnyango.csv")
competenceDat <- competenceDat[competenceDat$Ref %in% "Onyango 2020",]
bin <- binom.confint(x=competenceDat$NumDiss,n=competenceDat$NumInf,method="exact")
competenceDat$meanDiss <- bin$mean
competenceDat$lowerDiss <- bin$lower
competenceDat$upperDiss <- bin$upper
#****************************************************

#*****************parameter estimates from midgut infection*********************
source("3_INFfittingNLLParameterEstimates.R")
#*********************************************************************
#**************************Function to wrap fitting and simulations for iteration*************
fitToDataFunc <- function(x
                          ,pIFit
                          ,muVFit
                          ,prodRate=10
                          ,cellSpread=10^-3.5
                          ,escapeRate=0.05
                          ,dataSet){
  a<-x
  trace<-3
  #********initial parameter values****
  init.pars.fit <- c(   
    log_prodRate =   log(prodRate) 
    ,log_cellSpread = log(cellSpread)  
    ,log_escapeRate = log(escapeRate)   
  )
  #********optimise*******
  optim.vals <- optim(par = init.pars.fit
                      , objFXN
                      , fixed.params = virus_params(infRate=pIFit,muV=muVFit)
                      , dat = dataSet
                      , control = list(trace = trace, maxit = 200, reltol = 10^-7.5)
                      , method = "Nelder-Mead" # 
                      , hessian = T)
  
  return(list(optim.vals))
}
#********************************************************************************


aeg <- mclapply(1:100,fitToDataFunc
                ,dataSet=competenceDat[(cimpetenceDat$Moz %in% "Ae. aegypti"),]
                ,pIFit=aeginfRate
                ,muVFit=aegMuV
                ,prodRate=50
                ,cellSpread=10^-3.5
                ,escapeRate=0.09
)

saveRDS(aeg,"aegDissFits.rds")


alb <- mclapply(1:1,fitToDataFunc
                ,dataSet=competenceDat[(competenceDat$Moz %in% "Ae. albopictus"),]
                ,pIFit=albinfRate
                ,muVFit=aegMuV
                ,prodRate=10
                ,cellSpread=10^-6
                ,escapeRate=0.005)
saveRDS(alb,"albDissFits.rds")










#*******************************Simulate with fitted values**********************
# read in model fits fits
# read in model fits fits
alb <- readRDS("albDissFits220322.rds")
aeg <- readRDS("aegDissFits220322.rds")

# function to change to dataframe and then bind together
# function to change to dataframe and then bind together
listToDat <- function(output=alb, species="alb"){
  params <- lapply(1:length(output),function(x){
    temp <- output[[x]]
    temp <- temp[[1]][1]
    temp <- do.call(c,temp)
    names(temp) <- c("prodRate","cellSpread","escapeRate")
    return(exp(temp))
  })
  params <- do.call(rbind.data.frame,params)
  names(params) <- c("prodRate","cellSpread","escapeRate")
  params$species=species
  return(params)
}

alb <- listToDat()
aeg <- listToDat(output=aeg,species="aeg")

#********************************************
start_time <- Sys.time()

#length(albLI24Params[,1])

alb$run <- 1:length(alb[,1])

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)

albSims <- parRapply(cl,alb,function(y){
  source("5_DISSmodelFunc.R")
  source("5_dataForFitting.R")
  source("5_INFfittingNLLParameterEstimates.R")
  virus_params <- function( muV=aegMuV
                              ,infRate=albinfRate
                              ,prodRate =   as.numeric(y[1])
                              ,cellSpread = as.numeric(y[2])
                              ,escapeRate = as.numeric(y[3])   
                              ,cMax = 400            
                              ,hMax = 900             
  )
  return(as.list(environment()))
  
  sim.res <- sim.func(x=10^round(onyango$Conc.Min[1],0),parms=virus_params())
  
  run <- y[5]
  return(cbind.data.frame(sim.res,run))
  
})

stopCluster(cl)


end_time <- Sys.time()
timeTaken <- end_time - start_time
timeTaken

albSimsdf <- do.call(rbind,albSims)
albSimsdf$run <- as.factor(albSimsdf$run)

saveRDS(albSimsdf,"albSims.rds")

