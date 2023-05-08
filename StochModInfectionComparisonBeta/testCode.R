library(here)
library(binom)
library(emdbook)
library(EnvStats)
library(parallel)
library(bbmle)

#****model code****
source(here("StochModInfectionComparisonBeta//INFRepeatModelFunc.R"))
#***likelihood****
source(here("StochModInfectionComparison//INFfittingFuncsNLL.R"))

library(plyr)

virus_params <- function(   muV = 0.13
                            ,infRate = 10^-7.5
                            ,prodRate = 50             
                            ,cellSpread = 10^-4        
                            ,escapeRate = 0.0024        
                            ,cMax = 400                 
)
  return(as.list(environment()))


params <- virus_params()
# replicate experiments across virus concentrations 
cl <- makeCluster(detectCores()-1)
clusterEvalQ(cl, {library(adaptivetau)})
environment(repeatInfModel) <- .GlobalEnv
clusterExport(cl, varlist=c("infectionModel","repeatInfModel","params"),
              envir=environment())

test <- parLapply(cl,1:30,function(y){
  #set.seed(y)
  s <- repeatInfModel(x=params
                      ,startingVirus=c(10^2,10^3,10^4,10^5,10^6,10^6.5,10^7,10^8))
  return(s)
})
stopCluster(cl)


test <- do.call(rbind,test) 
# establish how many runs failed
fails <- length(test$num[test$num %in% NA])

stats <- lapply(unique(test$conc),function(z){
  temp <- test[test$conc %in% z,]
  indPs <- temp$num/temp$denom
  meanP <- sum(temp$num)/sum(temp$denom)
  varP <- sum(temp$denom)*(sum(temp$num)/sum(temp$denom))*(1-sum(temp$num)/sum(temp$denom))
  plot(1:length(temp$conc),indPs,ylim=c(0,1),ylab="Probability from 30 mosquito sims")
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
    return(c(temp$conc[1],betaEst,meanP,varP))
})

stats <- do.call(rbind.data.frame,stats)
names(stats) <- c("conc","shape1","shape2","mean","var")


# of course if all experimental replicates simulated result in either 0 or 1 there is no way to 
# fit a beta distribution to the resulting probabilities - for the first and last concentrations
# here there is certainty in either no infection or infectionacross all experiments when dose is low or high enough.
# so underlying probablility would be 0 or 1 with n
concs <- unique(test$conc)
z<-concs[4]
temp <- test[test$conc %in% z,]
indPs <- temp$num/temp$denom
meanP <- sum(temp$num)/sum(temp$denom)
plot(1:length(temp$conc),indPs,ylim=c(0,1))
ebeta(indPs) 

dbetas <- dbeta(x=indPs,shape1=0.2,shape2=10,log=T)

nllBeta <- function(logpar1,logpar2,dat){
  dbetas <- dbeta(x=dat,shape1=exp(logpar1),shape2=exp(logpar2),log=T)
  dbetas <- dbetas[!dbetas %in% Inf]
  dbetas <- dbetas[!dbetas<0]
  ll <- sum(dbetas)
  return(-ll)
}

# if try this in mle2 get non-finite finite difference value
# if exclude zeros:

test1 <- mle2(function(logpar1,logpar2){nllBeta(logpar1,logpar2,indPs)}             # fit first model
          ,start=list(logpar1=log(0.2),logpar2=log(10)))

# get quite different estimates


# for virus concentrations where there are experiments with zero but not all
# zeros perhaps these zeros are being omitted? the estimated distribution
# looks sensible do the zeros get omitted?
randomBeta <- rbeta(30,shape1=stats$shape1[4],shape2=stats$shape2[4])
plot(1:30,randomBeta,ylim=c(0,1))

randomBeta <- rbeta(30,shape1=exp(coef(test1)[1]),shape2=exp(coef(test1)[2]))
plot(1:30,randomBeta,ylim=c(0,1))

dbetabinom(0
           ,shape1=stats$shape1[4]
           ,shape2=stats$shape2[4]
           ,size=30
           ,log=F)
