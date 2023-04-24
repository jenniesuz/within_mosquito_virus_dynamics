library(here)
library(binom)
library(emdbook)
library(EnvStats)
library(parallel)
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

test <- parLapply(cl,1:100,function(y){
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
  plot(1:length(temp$conc),indPs,ylim=c(0,1))
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
names(stats) <- c("conc","shape1","shape2","mean")


# of course if all experimental replicates simulated result in either 0 or 1 there is no way to 
# fit a beta distribution to the resulting probabilities - for the first and last concentrations
# here there is certainty in either no infection or infection across all experiments when dose is low or high enough.
# so underlying probablility would be 0 or 1 with no variance - set both shape parameters to zero?



rbetabinom(100,shape1=stats$shape1[7],shape2=stats$shape2[7],size=30,prob=stats$mean[7])
rbetabinom(100,shape1=0,shape2=0,size=30)

dbetabinom(28,shape1=stats$shape1[7],shape2=stats$shape2[7],size=30)
dbinom(30,size=30,prob=0)



