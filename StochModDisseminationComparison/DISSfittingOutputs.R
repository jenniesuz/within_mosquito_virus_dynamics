library(plyr)
library(here)
library(parallel)
library(binom)
library(ggplot2)

source(here(".//StochModDisseminationComparison//DISSmodelFunc.R"))

#*************escape rate differs**************
output <- readRDS(here("rdsOutputs//DissModelFitSepAllParmsEscapeRateDiff220811.rds"))

-2*(-output$value) + 2*4  # 1178
# estimates
MLEfits <- output[[1]]
round(exp(MLEfits),6)
#CIs
fisherInfMatrix <- solve(-output$hessian) ## invert the Hessian, to estimate the covar-var matrix of parameter estimates
round(exp(MLEfits + 1.96*sqrt(diag(fisherInfMatrix))),6)
round(exp(MLEfits - 1.96*sqrt(diag(fisherInfMatrix))),6)



prodRate1 <- exp(MLEfits)[1]
prodRate2 <- exp(MLEfits)[1]
cellSpread1 <- exp(MLEfits)[2]
cellSpread2 <- exp(MLEfits)[2]
escapeRate1 <- exp(MLEfits)[3]
escapeRate2 <- exp(MLEfits)[4]


#******************data to fit to********************
competenceDat <- read.csv(here(".//Data//datCiotaOnyango.csv"))
competenceDat <- competenceDat[competenceDat$Ref %in% "Onyango 2020",]
bin <- binom.confint(x=competenceDat$NumDiss,n=competenceDat$NumInf,method="exact")
competenceDat$meanDiss <- bin$mean
competenceDat$lowerDiss <- bin$lower
competenceDat$upperDiss <- bin$upper
#****************************************************

#*************************************************************************************
#****run model with fitted values****

simAeg <- sim.func(x=10^competenceDat$Conc.Min[competenceDat$Moz %in% "Ae. aegypti"][1]
                   ,muV = 1.544027e-01
                   ,infRate = 2.338052e-09
                   ,prodRate = as.numeric(prodRate1)
                   ,cellSpread = as.numeric(cellSpread1)
                   ,escapeRate = as.numeric(escapeRate1)
                   ,cMax = 400 
                   ,hMax = 900)

#**********first remove simulations where an infection didn't occur**********************
inf <- ddply(simAeg,.(run),summarise,occurred=sum(inf))  # establish if infection occurred
inf$occurred[inf$occurred>0] <- 1
noInf <- inf$run[inf$occurred == 0]                         
simAeg <- simAeg[!simAeg$run %in% noInf,]
#*******************************************************************************

simAeg2 <- dissSummaryFunc2(simAeg)


simAlb <- sim.func(x=10^competenceDat$Conc.Min[competenceDat$Moz %in% "Ae. albopictus"][1]
                   ,muV = 1.544027e-01
                   ,infRate = 5.494484e-08 
                   ,prodRate =   as.numeric(prodRate2)
                   ,cellSpread = as.numeric(cellSpread2)
                   ,escapeRate =  as.numeric(escapeRate2)
                   ,cMax = 400 
                   ,hMax = 900)

#**********first remove simulations where an infection didn't occur**********************
inf <- ddply(simAlb,.(run),summarise,occurred=sum(inf))  # establish if infection occurred
inf$occurred[inf$occurred>0] <- 1
noInf <- inf$run[inf$occurred == 0]                         
simAlb <- simAlb[!simAlb$run %in% noInf,]
#*******************************************************************************

simAlb2 <- dissSummaryFunc2(simAlb)

#*************plot*************
names(simAeg2) <- c("DPIDissorTrans","meanDiss")
names(simAlb2) <- c("DPIDissorTrans","meanDiss")

ggplot(competenceDat) +
  geom_line(data=simAeg2,aes(x=DPIDissorTrans,y=meanDiss),linetype=3,alpha=0.7,col="royalblue4") +
  geom_line(data=simAlb2,aes(x=DPIDissorTrans,y=meanDiss),linetype=3,alpha=0.7,col="dodgerblue") +
  geom_errorbar(aes(x=DPIDissorTrans,ymin=lowerDiss,ymax=upperDiss),alpha=0.5) +
  geom_point(aes(x=DPIDissorTrans,y=meanDiss,fill=factor(Moz)),shape=21) +
  scale_fill_manual(values=c("royalblue4","dodgerblue")) +
  labs(fill="",title="B") +
  xlab("Days post blood meal") +
  ylab("Proportion of infected mosquitoes \n with a disseminated infection") +
  theme_set(theme_bw())  +    
  theme(panel.border = element_blank()                   
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=8)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.5,"line")   # here you can amend legend size and position
        ,legend.background = element_blank()
        ,legend.text=element_text(size=6)
        ,strip.text=element_text(size=5)
        ,legend.position="none"
        ,strip.background = element_rect(fill="white",color="white")
  )  
#************************************













#*************all params**************
output <- readRDS(here("rdsOutputs//DissModelFitSepAllParms220804.rds"))

-2*(-output$value) + 2*6  # 1062
MLEfits <- output[[1]]
round(exp(MLEfits),6)
fisherInfMatrix <- solve(-output$hessian) ## invert the Hessian, to estimate the covar-var matrix of parameter estimates
exp(MLEfits + 1.96*sqrt(diag(fisherInfMatrix)))
exp(MLEfits - 1.96*sqrt(diag(fisherInfMatrix)))


# estimates
prodRate1 <- exp(MLEfits)[1]
prodRate2 <- exp(MLEfits)[2]
cellSpread1 <- exp(MLEfits)[3]
cellSpread2 <- exp(MLEfits)[4]
escapeRate1 <- exp(MLEfits)[5]
escapeRate2 <- exp(MLEfits)[6]


#******************data to fit to********************
competenceDat <- read.csv(here(".//Data//datCiotaOnyango.csv"))
competenceDat <- competenceDat[competenceDat$Ref %in% "Onyango 2020",]
bin <- binom.confint(x=competenceDat$NumDiss,n=competenceDat$NumInf,method="exact")
competenceDat$meanDiss <- bin$mean
competenceDat$lowerDiss <- bin$lower
competenceDat$upperDiss <- bin$upper
#****************************************************

#*************************************************************************************
#****run model with fitted values****

simAeg <- sim.func(x=10^competenceDat$Conc.Min[competenceDat$Moz %in% "Ae. aegypti"][1]
                   ,muV = 1.544027e-01
                   ,infRate = 2.338052e-09
                   ,prodRate =  100# as.numeric(prodRate1)
                   ,cellSpread = 10^-3.5#as.numeric(cellSpread1)
                   ,escapeRate = 0.005#as.numeric(escapeRate1)
                   ,cMax = 400 
                   ,hMax = 900)

#**********first remove simulations where an infection didn't occur**********************
inf <- ddply(simAeg,.(run),summarise,occurred=sum(inf))  # establish if infection occurred
inf$occurred[inf$occurred>0] <- 1
noInf <- inf$run[inf$occurred == 0]                         
simAeg <- simAeg[!simAeg$run %in% noInf,]
#*******************************************************************************

simAeg2 <- dissSummaryFunc2(simAeg)


simAlb <- sim.func(x=10^competenceDat$Conc.Min[competenceDat$Moz %in% "Ae. albopictus"][1]
                   ,muV = 1.544027e-01
                   ,infRate = 5.494484e-08 
                   ,prodRate =   100#as.numeric(prodRate2)
                   ,cellSpread = 10^-3.5#as.numeric(cellSpread2)
                   ,escapeRate =  0.0001 #as.numeric(escapeRate2)
                   ,cMax = 400 
                   ,hMax = 900)

#**********first remove simulations where an infection didn't occur**********************
inf <- ddply(simAlb,.(run),summarise,occurred=sum(inf))  # establish if infection occurred
inf$occurred[inf$occurred>0] <- 1
noInf <- inf$run[inf$occurred == 0]                         
simAlb <- simAlb[!simAlb$run %in% noInf,]
#*******************************************************************************

simAlb2 <- dissSummaryFunc2(simAlb)

#*************plot*************
names(simAeg2) <- c("DPIDissorTrans","meanDiss")
names(simAlb2) <- c("DPIDissorTrans","meanDiss")

ggplot(competenceDat) +
  geom_line(data=simAeg2,aes(x=DPIDissorTrans,y=meanDiss),linetype=3,alpha=0.7,col="royalblue4") +
  geom_line(data=simAlb2,aes(x=DPIDissorTrans,y=meanDiss),linetype=3,alpha=0.7,col="dodgerblue") +
  geom_errorbar(aes(x=DPIDissorTrans,ymin=lowerDiss,ymax=upperDiss),alpha=0.5) +
  geom_point(aes(x=DPIDissorTrans,y=meanDiss,fill=factor(Moz)),shape=21) +
  scale_fill_manual(values=c("royalblue4","dodgerblue")) +
  labs(fill="",title="B") +
  xlab("Days post blood meal") +
  ylab("Proportion of infected mosquitoes \n with a disseminated infection") +
  theme_set(theme_bw())  +    
  theme(panel.border = element_blank()                   
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=8)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.5,"line")   # here you can amend legend size and position
        ,legend.background = element_blank()
        ,legend.text=element_text(size=6)
        ,strip.text=element_text(size=5)
        ,legend.position="none"
        ,strip.background = element_rect(fill="white",color="white")
  )  
#************************************