library(EasyABC)
library(binom)
library(parallel)


#******************data to fit to********************
competenceDat <- read.csv(here::here(".//Data//datCiotaOnyango.csv"))
competenceDat <- competenceDat[competenceDat$Ref %in% "Ciota 2017",]
bin <- binom.confint(x=competenceDat$NumInf,n=competenceDat$ITotal,method="exact")
competenceDat$meanInf <- bin$mean
competenceDat$lowerInf <- bin$lower
competenceDat$upperInf <- bin$upper
#****************************************************
aeg <- competenceDat[competenceDat$Moz %in% "Ae. aegypti",]
modAeg <- glm(cbind(NumInf,ITotal-NumInf) ~ ConcMax,weights=ITotal,family="binomial",data=aeg)
summary(modAeg)
coef(modAeg)

# if we assume these are the 'true' parameter values and we want to estimate them

toyMod <- function(x=c(par1=-8,par2=1)
                   ,virusConcs=competenceDat$ConcMax[competenceDat$Moz %in% "Ae. aegypti"]
){
  y <- tryCatch( {locfit::expit(x[1] + x[2]*virusConcs)} , error = function(e) {
    print(e)
    return(c(NA,NA))
  })
  if(is.na(y[1])==F){ 
    num1 <- rbinom(n=30,prob=y[1],size=1)  # simulate data based on above
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


#**************************************************
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

plot(ABC_Beaumont$param[,1],ABC_Beaumont$weight)
plot(ABC_Beaumont$param[,2],ABC_Beaumont$weight)



startTime <- Sys.time()
ABC_Lenormand <- ABC_sequential(method="Lenormand"
                               , model=toyMod
                               , prior=priorParams
                               , nb_simul=1000
                               , summary_stat_target=sum_stat_obs
                               ,verbose=F
                               ,use_seed=F
                               ,alpha=0.9
                               ,prior_test="X2>X1")

endTime <- Sys.time()
endTime - startTime
# 1.6 mins when alpha = 0.5
# 

hist(ABC_Lenormand$param[,1])
hist(ABC_Lenormand$param[,2])

plot(ABC_Lenormand$param[,1],ABC_Lenormand$weight)
plot(ABC_Lenormand$param[,2],ABC_Lenormand$weight)

