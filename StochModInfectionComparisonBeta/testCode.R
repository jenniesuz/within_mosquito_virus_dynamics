library(here)
library(binom)
library(emdbook)
#****model code****
source(here("StochModInfectionComparisonBeta//INFRepeatModelFunc.R"))
#***likelihood****
source(here("StochModInfectionComparison//INFfittingFuncsNLL.R"))

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
  totI <- sum(temp$num,na.rm=T)
  totD <- sum(temp$denom,na.rm=T)
  meanProp <- totI/totD
  meanCount <- mean(temp$num,na.rm=T)
  indPs <- temp$num/temp$denom
  varProp <- sum((indPs - meanProp)^2)/(length(indPs - 1))
  varCount <- sum((temp$num - meanCount)^2)/length(temp$num - 1)
  plot(1:length(temp$conc),indPs,ylim=c(0,1))
    return(c(temp$conc[1],meanProp,meanCount,varProp,varCount,totI))
})

stats <- do.call(rbind.data.frame,stats)
names(stats) <- c("conc","meanProp","meanCount","varProp","varCount","totI")



test6 <- rbinom(100,30,stats$meanProp[6])
plot(1:100,test6/30,ylim=c(0,1))

a <- function(mean,var){
  return( ( ((1-mean)/var^2) - 1/mean )*mean^2 )
}

b <- function(a,mean){
  return(  a*(1/mean - 1 ))
}

dbetabinom()

rbetabinom(30,stats$meanProp[3],size=100,)