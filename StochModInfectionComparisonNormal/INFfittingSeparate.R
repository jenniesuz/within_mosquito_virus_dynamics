library(here)
library(binom)
library(ggplot2)
library(bbmle)
library(parallel)
library(adaptivetau)
library(mvtnorm)
#****model code****
source(here("StochModInfectionComparisonNormal//INFmodelFunc.R"))

source(here("StochModInfectionComparisonNormal//INFRepeatModelFunc.R"))
#***likelihood****
source(here("StochModInfectionComparisonNormal//INFfittingFuncsNLL.R"))

#******************data to fit to********************
competenceDat <- read.csv(here("StochModInfectionComparisonNormal//datCiotaOnyango.csv"))
competenceDat <- competenceDat[competenceDat$Ref %in% "Ciota 2017",]
bin <- binom.confint(x=competenceDat$NumInf,n=competenceDat$ITotal,method="exact")
competenceDat$meanInf <- bin$mean
competenceDat$lowerInf <- bin$lower
competenceDat$upperInf <- bin$upper
#****************************************************

# test for difference between species in dose-response curve using logistic regression
noSppMod <- glm(NumInf/ITotal ~ ConcMax,data=competenceDat,family="binomial",weights=ITotal)
SppMod <- glm(NumInf/ITotal ~ ConcMax+Moz,data=competenceDat,family="binomial",weights=ITotal)
SppModInt <- glm(NumInf/ITotal ~ ConcMax*Moz,data=competenceDat,family="binomial",weights=ITotal)
AIC(noSppMod)
AIC(SppMod)
AIC(SppModInt)



#**********************get approx starting parameters*********************
ggplot(competenceDat) +
  geom_point(aes(x=ConcMax,y=meanInf,col=Moz))

virus_params <- function(   muV = 0.1
                            ,infRate1 = 10^-8
                            ,prodRate1 = 10            # note these parameters don't matter for infection
                            ,cellSpread1 = 10^-4        
                            ,escapeRate1 = 0.05         #
                            ,cMax = 400                 #
)
  return(as.list(environment()))

testAeg <- repeatInfModel(x=virus_params(infRate1=10^-8)
 ,startingVirus=10^competenceDat$ConcMax[competenceDat$Moz %in% "Ae. aegypti"]
)
testAeg

#************************parameters for two 'treatments'****************************
virus_params <- function(   muV = 0.1
                            ,infRate1 = 10^-8
                            ,infRate2 = 10^-8
                            ,prodRate1 = 1              # note these parameters don't matter for infection
                            ,cellSpread1 = 10^-4        
                            ,escapeRate1 = 0.05         #
                            ,cMax = 400                 #
)
  return(as.list(environment()))
#*****************************************************************************
init.pars.fit <- c(
  log_infRate1=log(10^-9)
  ,log_infRate2=log(10^-7.5)
)

objFXN(fit.params=init.pars.fit                                                          ## paramters to fit
                   , fixed.params =virus_params()                                      ## fixed paramters
                   , dat=competenceDat
                   , nSimulations=30)

#************curve******************
start <- Sys.time()
testLik30 <- lapply(10^-seq(7,9,0.1)
                    ,function(x){
  loglik <- nll.binom(parms=virus_params(muV=0.1
                                         ,infRate1=x
                                         ,infRate2=x)
                      ,dat=competenceDat
                      ,nSims=30) 
  return(loglik)
})

testLik30 <- do.call(c,testLik30)
end <- Sys.time()
end-start

dat <- cbind.data.frame(beta=10^-seq(7,9,0.1),ll=testLik30)

ggplot(dat) +
  geom_line(aes(x=log10(beta),y=ll)) +
  ylim(0,25)



#******************************************

start <- Sys.time()
  trace<-3
  #********initial parameter values****
  init.pars.fit <- c(
    log_infRate1= -20.48997
    ,log_infRate2=-19.00023  
  )
  
  #********optimise*******
  optim.vals <- optim(par = init.pars.fit
                      , objFXN
                      , fixed.params = virus_params()
                      , dat = competenceDat
                      , nSimulations = 30
                      , control = list(trace = trace
                                       ,abstol=0.05
                                       ,reltol=0.05
                                       #,maxit=200
                                      )
                      , method ="Nelder-Mead" #"SANN" #
                     )
  
  end <- Sys.time()
  end-start
  # save fitted parameters
  #AIC
  -2*(-optim.vals$value) + 2*2
  
  #saveRDS(optim.vals,"INFModelFitDiffParms110523NM.rds")
  
# same parms started 10^-7
# 1st SANN -16.85324  
# 2nd SANN -16.67233 
# NM rel tol 0.05 converged -16.67233 $value [1] 10.85173   #23.703


# diff parms started 10^-7
# 1st SANN  log_infRate1 log_infRate2 
 # -19.19664    -16.90687
# 2nd SANN
 # log_infRate1 log_infRate2 
#  -20.48997    -19.00023 
# NM rel to 0.05 converged log_infRate1 log_infRate2 
                           # -20.64605    -19.00823    # 15.964


#****************************
nSims <- 100
parmsAeg <- c(0.1
              ,exp(-19.19664)
              ,virus_params()$prodRate           
              ,virus_params()$cellSpread        
              ,virus_params()$escapeRate          
              ,virus_params()$cMax)
datAeg <- competenceDat[competenceDat$Moz %in% "Ae. aegypti",]

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

linesAeg <- 