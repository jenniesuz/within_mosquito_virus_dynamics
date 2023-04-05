library(EasyABC)
library(binom)
source(here::here("StochModInfectionComparisonABCPar//INFRepeatModelFunc.R"))


#******************data to fit to********************
competenceDat <- read.csv(here::here(".//Data//datCiotaOnyango.csv"))
competenceDat <- competenceDat[competenceDat$Ref %in% "Ciota 2017",]
bin <- binom.confint(x=competenceDat$NumInf,n=competenceDat$ITotal,method="exact")
competenceDat$meanInf <- bin$mean
competenceDat$lowerInf <- bin$lower
competenceDat$upperInf <- bin$upper
#****************************************************

modAeg <- glm(matrix(c(NumInf,ITotal-NumInf),ncol=2) ~ ConcMax
              ,family="binomial"
              ,data=competenceDat[competenceDat$Moz %in% "Ae. aegypti",])
coef(modAeg)


toyMod <- function(x=c(par1=-8,par2=1)
          ,virusConcs=competenceDat$ConcMax[competenceDat$Moz %in% "Ae. aegypti"]
  ){
  y <- tryCatch( {locfit::expit(x[1] + x[2]*virusConcs)} , error = function(e) {
    print(e)
    return(c(NA,NA))
  })
  if(is.na(y[1])==F){ 
  num1 <- rbinom(n=30,prob=y[1],size=1)
  num2 <- rbinom(n=30,prob=y[2],size=1)
  num3 <- rbinom(n=30,prob=y[3],size=1)
  num4 <- rbinom(n=30,prob=y[4],size=1)
  simDat <- cbind.data.frame(event=c(num1,num2,num3,num4)
                             ,conc=c(rep(virusConcs[1],30)
                                     ,rep(virusConcs[2],30)
                                     ,rep(virusConcs[3],30)
                                     ,rep(virusConcs[4],30))
  )
  modGLM <- glm(event ~ conc
                ,family="binomial"
                ,data=simDat)
  return(as.numeric(coef(modGLM)))
  }else{ return(c(-999,-999))}
}

priorParams <- list(c("unif",-20,-2)
                    ,c("unif",0,20))
sum_stat_obs <- as.numeric(coef(modAeg))
tolerance <- c(2,1)

startTime <- Sys.time()

ABC_Beaumont <- ABC_sequential(method="Beaumont"
                             , model=toyMod
                             , prior=priorParams
                             , nb_simul=1000
                             , summary_stat_target=sum_stat_obs
                             , tolerance_tab=tolerance
                             ,verbose=F
                             ,use_seed=F
                             ,prior_test="X2>X1")

endTime <- Sys.time()
endTime - startTime

hist(ABC_Beaumont$param[,1])
hist(ABC_Beaumont$param[,2])

hist(ABC_Beaumont$stats[,1])
hist(ABC_Beaumont$stats[,2])

hist(ABC_Beaumont$weights)


plot(ABC_Beaumont$param[,1],ABC_Beaumont$weight)
plot(ABC_Beaumont$param[,2],ABC_Beaumont$weight)

