library(gridExtra)
library(ggplot2)
library(plyr)
library(ggfan)
library(binom)
source("4_DISSmodelFunc.R")
#***************Midgut Simulations**********************

aegInfSims <- readRDS("aegInfSims.rds")
albInfSims  <- readRDS("albInfSims.rds")

aegInfSims <- do.call(rbind.data.frame,aegInfSims)
aegInfSims$Moz <- "Ae. aegypti"
aegInfSims$parmsComb <- paste(aegInfSims$probInf,aegInfSims$muV)
albInfSims <- do.call(rbind.data.frame,albInfSims)
albInfSims$Moz <- "Ae. albopictus"
albInfSims$parmsComb <- paste(albInfSims$probInf,albInfSims$muV)

sims <- rbind.data.frame(aegInfSims,albInfSims)

#******************data to fit to********************
dat<- read.csv("datCiotaOnyango.csv")
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
#*
infFitsPlot <- ggplot(dat[dat$Ref %in% "Ciota 2017",]) +
  geom_line(data=sims,aes(x=log10(Conc.Min),y=meanInf,col=Moz,group=parmsComb),linetype=3,alpha=0.8) +
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




#***************Diss Simulations**********************


aegDissSims <- read.csv("aegDissSummary220418.csv")
names(aegDissSims)[2] <- "DPIDissorTrans"
names(aegDissSims)[3] <- "meanDiss"


albDissSims <- read.csv("albDissSummary220418.csv")
names(albDissSims)[2] <- "DPIDissorTrans"
names(albDissSims)[3] <- "meanDiss"

aegDissSims$Moz <- ("Ae. aegypti")
albDissSims$Moz <-  ("Ae. albopictus")

sims <- rbind.data.frame(aegDissSims,albDissSims) 

dissFitsPlot <- ggplot(dat[dat$Ref %in% "Onyango 2020",]) +
  geom_line(data=sims[sims$Moz %in% "Ae. aegypti",],aes(x=DPIDissorTrans,y=meanDiss,group=mainRun),linetype=3,alpha=0.7,col="royalblue4") +
  geom_line(data=sims[sims$Moz %in% "Ae. albopictus",],aes(x=DPIDissorTrans,y=meanDiss,group=mainRun),linetype=3,alpha=0.7,col="dodgerblue") +
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
  

pdf(file="fig_midgutHemocoelModelFits.pdf",width=6,height=4)
grid.arrange(infFitsPlot,dissFitsPlot,ncol=2)
dev.off()