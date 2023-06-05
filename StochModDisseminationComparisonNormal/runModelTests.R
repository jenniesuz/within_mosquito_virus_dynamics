
library(here)
source(here(".//StochModDisseminationComparisonNormal//InfDissRepeatModelFunc.R"))
source(here(".//stochasticModels//InfDissModelFuncs.R"))


library(grid)
library(gridExtra)
library(ggplot2)

#******************data to fit to********************
competenceDat <- read.csv(here(".//Data//datCiotaOnyango.csv"))
competenceDat <- competenceDat[competenceDat$Ref %in% "Onyango 2020",]
bin <- binom.confint(x=competenceDat$NumDiss,n=competenceDat$NumInf,method="exact")
competenceDat$meanDiss <- bin$mean
competenceDat$lowerDiss <- bin$lower
competenceDat$upperDiss <- bin$upper
#****************************************************

virus_params <- function(   muV = 0.1
                            ,infRate = 10^-7 #10^-8.8
                            ,prodRate = 25          # note these parameters don't matter for infection
                            ,cellSpread = 0.00005       
                            ,escapeRate = 0.005  #
                            ,cMax = 400  
                            ,hMax = 900             #
)
  return(as.list(environment()))

#3.673309e-09
#***************test code**************
startTime <- Sys.time()
test <- repeatModel(virus_params(infRate=10^-8),startingVirus=10^9)
endTime <- Sys.time()
endTime - startTime

test <- do.call(rbind,test)

ggplot(test) +
  geom_line(aes(x=time,y=Hv,group=run,col=run))

#* calculate proportions 
finalTime <- 400

maxTime <-  max(test$time)
if(maxTime < finalTime){
  restTimes <- seq((maxTime+0.3),finalTime,0.3)
  
  zerosDat <- cbind.data.frame(restTimes
                               ,rep(0,length(restTimes))
                               ,rep(0,length(restTimes))
                               ,rep(0,length(restTimes))
                               ,rep(0,length(restTimes))
                               ,rep(0,length(restTimes)))
  names(zerosDat) <- c("t","Gv","Mci","Mv","Hv","Hci")
  out <- rbind.data.frame(out,zerosDat)
}


#************Midgut infection***********
sum(test$inf)/length(test$inf)
 
#*******************************

testProps <- dissSummaryFunc(test)

ggplot(testProps) +
  geom_point(data=competenceDat,aes(x=DPIDissorTrans,y=meanDiss)) +
  geom_errorbar(data=competenceDat,aes(x=DPIDissorTrans,ymin=lowerDiss,ymax=upperDiss)) +
  geom_line(aes(x=time,y=proportionDisseminated)) +
  xlim(0,14)


