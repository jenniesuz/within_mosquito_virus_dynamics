library(here)
library(binom)
library(parallel)
library(bbmle)
library(mvtnorm)

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
  s$run <- y
  return(s)
})
stopCluster(cl)

test <- do.call(rbind,test) 
# establish how many runs failed
fails <- length(test$num[test$num %in% NA])

stats <- lapply(unique(test$run),function(z){
  temp <- test[test$run %in% z,]
  mod <- glm(temp$num/temp$denom~temp$conc,family="binomial",weights=temp$denom)
  return(c(temp$run[1],as.vector(coef(mod))))
})

stats <- do.call(rbind.data.frame,stats)
names(stats) <- c("run","p1","p2")


#****************************************************
cl <- makeCluster(detectCores()-1)
clusterEvalQ(cl, {library(adaptivetau)})
environment(repeatInfModel) <- .GlobalEnv
clusterExport(cl, varlist=c("infectionModel","repeatInfModel","params"),
              envir=environment())

test <- parLapply(cl,1:30,function(y){
  #set.seed(y)
  s <- repeatInfModel(x=params
                      ,startingVirus=c(10^2,10^3,10^4,10^5,10^6,10^6.5,10^7,10^8))
  s$run <- y
  return(s)
})
stopCluster(cl)

test <- do.call(rbind,test) 
# establish how many runs failed
fails <- length(test$num[test$num %in% NA])

stats2 <- lapply(unique(test$run),function(z){
  temp <- test[test$run %in% z,]
  mod <- glm(temp$num/temp$denom~temp$conc,family="binomial",weights=temp$denom)
  return(c(temp$run[1],as.vector(coef(mod))))
})

stats2 <- do.call(rbind.data.frame,stats2)
names(stats2) <- c("run","p1","p2")



dmvnorm(stats2[1,2:3]
        ,mean=c(mean(stats[,2])
                ,mean(stats[,3]))
        ,sigma=cov(stats[,2:3]),log=T)
