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

priorParams <- list(c("unif",0.01,0.04)
                   ,c("unif",10^-10,10^-6))
sum_stat_obs <- as.numeric(coef(modAeg))
tolerance <- c(2,1)

startTime <- Sys.time()

#cl <- makeCluster(detectCores()-1)
#clusterEvalQ(cl, {library(adaptivetau)})
#environment(infectionModel) <- .GlobalEnv
#clusterExport(cl, varlist=c("infectionModel","virusConcs","x"),
#              envir=environment())


ABC_Beaumont <- ABC_sequential(method="Beaumont"
                             , model=repeatInfModel
                             , prior=priorParams
                             , nb_simul=10
                             , summary_stat_target=sum_stat_obs
                             , tolerance_tab=tolerance
                             ,verbose=F
                             ,use_seed=T
                             ,n_cluster=detectCores()-1)

endTime <- Sys.time()
endTime - startTime

hist(ABC_Beaumont$param[,1])
hist(ABC_Beaumont$param[,2])
hist(ABC_Beaumont$stats)
hist(ABC_Beaumont$weights)
plot(ABC_Beaumont$param[,2],ABC_Beaumont$weight)

