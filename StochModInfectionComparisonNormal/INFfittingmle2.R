library(bbmle)
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



dat=competenceDat[competenceDat$Moz %in% "Ae. albopictus",]


seqs <- seq(6,10,0.05)
paramSeq <- 10^-seqs
paramSeq <- rep(paramSeq,30)
paramSeq <- sort(paramSeq)

start <- Sys.time()

lls <- lapply(paramSeq,function(x){
  ll <- nll.binom.mle2(loginfRate=log(x),logmuV=log(0.1),dat=competenceDat[competenceDat$Moz %in% "Ae. albopictus",]
)
  return(ll)
})

end <- Sys.time()
end-start


llsd <- unlist(lls)
llsdp <- cbind.data.frame(llsd,paramSeq)

saveRDS(llsdp,"albLLRange.230603")

meanllsdp <- ddply(llsdp,.(paramSeq),summarise,mean(llsd,na.rm=T))
names(meanllsdp) <- c("paramSeq","llsd")
plot(log10(meanllsdp$paramSeq),-meanllsdp$llsd,xlim=c(-10,-6),ylim=c(-20,-6),type="l",col="red")
abline(h=max(-meanllsdp$llsd,na.rm=T)-1.92)
#************************************************************

# Ae. aeg c. 10^-8.5 and 10^-9 then c. 10^-9.4 to 10^-9.8
# min appears c. 10^-8.9


# Ae. alb c. 10^-6.8 (v. small) and 10^-8 to c. 10^-9.4

virus_params <- function(   muV = 0.1
                            ,infRate1 = 10^-8
                            ,prodRate1 = 1           # note these parameters don't matter for infection
                            ,cellSpread1 = 10^-4        
                            ,escapeRate1 = 0.05         #
                            ,cMax = 400                 #
)
  return(as.list(environment()))




testAeg <- repeatInfModel(x=virus_params(infRate1=10^-8)
                          ,startingVirus=10^competenceDat$ConcMax[competenceDat$Moz %in% "Ae. aegypti"]
)

newData <- cbind.data.frame("conc"=seq(1,10,0.1))

modAeg <- glm(num/denom~conc,family="binomial",weights=denom,data=testAeg)
predAeg<-predict(modAeg,type="response",newdata=newData)
predAegd <- cbind.data.frame("conc"=seq(1,10,0.1),"predVals"=predAeg)


ggplot(dat) +
  geom_point(aes(x=ConcMax,y=meanInf)) +
  geom_errorbar(aes(x=ConcMax,ymin=lowerInf,ymax=upperInf)) +
  geom_line(data=predAegd,aes(x=conc,y=predVals)) +
  xlim(3,10)




newData <- cbind.data.frame("ConcMax"=seq(1,10,0.1))

modAegData <- glm(NumInf/ITotal~ConcMax,family="binomial",weights=ITotal,data=dat)
predAegData <-predict(modAegData,type="response",newdata=newData)
predAegdData <- cbind.data.frame("conc"=seq(1,10,0.1),"predVals"=predAegData)




ggplot(dat) +
  geom_point(aes(x=ConcMax,y=meanInf)) +
  geom_errorbar(aes(x=ConcMax,ymin=lowerInf,ymax=upperInf)) +
  geom_line(data=predAegdData,aes(x=conc,y=predVals)) +
  xlim(3,10)


