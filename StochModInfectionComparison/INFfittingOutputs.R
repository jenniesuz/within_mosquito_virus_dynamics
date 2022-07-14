
library(here)
library(parallel)
library(binom)
library(ggplot2)
library(GillespieSSA)

source(here(".//StochModInfectionComparison//INFmodelFunc.R"))
source(here(".//StochModInfectionComparison//INFfittingFuncsNLL.R"))


output <- readRDS(here("rdsOutputs//INFModelFitSepParms220712.rds"))
# estimates
infParms <- exp(output[[1]])

log10(infParms)

-2*(-output$value) + 2*3 # 57.9





#***************************************
#******************data to fit to********************
dat<- read.csv(here("Data\\datCiotaOnyango.csv"))
#****************************************************
bin1 <- binom.confint(x=dat$NumInf,n=dat$ITotal,method="exact")
bin2 <- binom.confint(x=dat$NumDiss,n=dat$NumInf,method="exact")
dat$meanInf <- bin1$mean
dat$lowerInf <- bin1$lower
dat$upperInf <- bin1$upper
dat$meanDiss <- bin2$mean
dat$lowerDiss <- bin2$lower
dat$upperDiss <- bin2$upper
#***************************************************
library(plyr)
#*************************************************************************************
#****run model with fitted values****
doses <- c(10^3,10^3.5,10^4,10^4.5,10^5,10^5.5,10^6,10^6.5,10^7,10^7.5,10^8,10^8.5,10^9,10^9.5,10^10)


simAeg <- mclapply(doses,
                 sim.func
                 ,muV = as.numeric(infParms[1])
                 ,infRate = as.numeric(infParms[2]))

doseSimAeg <- do.call(rbind.data.frame,simAeg)
infDatAeg <- modOutFunc(doseSimAeg)

simAlb <- mclapply(doses,
                   sim.func
                   ,muV = as.numeric(infParms[1])
                   ,infRate = as.numeric(infParms[3]))

doseSimAlb <- do.call(rbind.data.frame,simAlb)
infDatAlb <- modOutFunc(doseSimAlb)

infDatAeg$Moz <- "Ae. aegypti"
infDatAlb$Moz <- "Ae. albopictus"

infDat <- rbind.data.frame(infDatAeg,infDatAlb)

ggplot(dat[dat$Ref %in% "Ciota 2017",]) +
  geom_line(data=infDat,aes(x=log10(Conc.Min),y=meanInf,col=Moz),linetype=3,alpha=0.8) +
  geom_errorbar(aes(x=Conc.Min,ymin=lowerInf,ymax=upperInf),alpha=0.5) +
  geom_point(aes(x=Conc.Min,y=meanInf,fill=Moz),shape=21)+
  # geom_point(data=sims,aes(x=log10(Conc.Min),y=meanInf,col=Moz),size=2)  +
  xlab(expression(paste("Virus concentration (log"[10]," PFU/ml)"))) +
  labs(fill="",title="A") +
  scale_fill_manual(values=c("royalblue4","dodgerblue")) +
  scale_colour_manual(values=c("royalblue4","dodgerblue"),guide="none") +
  ylab("Proportion of mosquitoes with a midgut infection") +
  theme_set(theme_bw())  +    
  theme(panel.border = element_blank()                   
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=8)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.5,"line")  
        ,legend.background = element_blank()
        ,legend.text=element_text(size=6)
        ,strip.text=element_text(size=5)
        ,legend.position=c(0.2,0.8)
  )   

