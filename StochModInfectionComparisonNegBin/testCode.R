library(here)
library(binom)
library(parallel)
library(poisbinom)
#****model code****
source(here("StochModInfectionComparisonNegBin//INFRepeatModelFunc.R"))
#***likelihood****
#source(here("StochModInfectionComparison//INFfittingFuncsNLL.R"))

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


## Binomial probabilities
#pp <- runif(30)
## PMF
#dpoisbinom(c(2), pp)

stats <- lapply(unique(test$conc),function(z){
  temp <- test[test$conc %in% z,]
  meanN <- round(mean(temp$num),0)
  varN <- var(temp$num)
  sizeN <- meanN^2/(varN-meanN)
    return(list(temp$conc[1],meanN,varN,sizeN))
})

stats <- do.call(rbind.data.frame,stats)
names(stats) <- c("conc","mean","var")


dnorm(0,mean=stats[1,2],sd=sqrt(stats[1,3]),log=T)
