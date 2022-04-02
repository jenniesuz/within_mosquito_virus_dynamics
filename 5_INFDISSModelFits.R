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
bin2 <- binom.confint(x=dat$NumDiss,n=dat$DTotal,method="exact")
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
  scale_colour_manual(values=c("royalblue4","dodgerblue")) +
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
        ,legend.position="none"
  )   




#***************Diss Simulations**********************
aegDissSims <- readRDS("aegDissSims.rds")
albDissSims <- readRDS("albDissSims.rds")

names(aegDissSims)[7] <- "subRun"
names(albDissSims)[7] <- "subRun"

#**********first remove simulations where an infection didn't occur**********************
inf <- ddply(aegDissSims,.(run,conc),summarise,occurred=sum(inf))  # establish if infection occurred
inf$occurred[inf$occurred>0] <- 1
inf$runConc <- paste(inf$run,inf$conc)
noInf <- inf$runConc[inf$occurred == 0]                         # note run and concs where infection didn't occur

aegDissSims$runConc <- paste(aegDissSims$run,aegDissSims$conc)           # remove these from the dissemination dataset
aegDissSims <- aegDissSims[!aegDissSims$runConc %in% noInf,]
#*******************************************************************************

aegDissSims2 <- lapply(unique(aegDissSims$run),function(y){
  temp <- aegDissSims[aegDissSims$run %in% y,]
    diss <- dissSummaryFunc(temp)
  diss$run <- y
return(diss)
})
aegDissSims2 <- do.call(rbind.data.frame,aegDissSims2)

write.csv(aegDissSims2,"aegDissSims.csv")


aegDissSims <- read.csv("aegDissSims.csv")
names(aegDissSims)[3] <- "DPIDissorTrans"
names(aegDissSims)[2] <- "meanDiss"

#********************
albDissSims2 <- lapply(unique(albDissSims$run),function(y){
  temp <- albDissSims[albDissSims$run %in% y,]
  diss <- dissSummaryFunc(temp)
  diss$run <- y
  return(diss)
})
albDissSims2 <- do.call(rbind.data.frame,albDissSims2)

write.csv(albDissSims2,"albDissSims.csv")


#**********

alb <- read.csv("albDissSims.csv")
aeg<- read.csv("aegDissSims.csv")

alb$species <- ("Ae. albopictus")
aeg$species <- ("Ae. aegypti")

sims <- rbind.data.frame(alb,aeg) 


#****************************************************
names(sims)[2] <- "DPIDissorTrans"
names(sims)[3] <- "meanDiss"
names(sims)[5] <- "Moz"


dissFitsPlot <- ggplot(dat[dat$Ref %in% "Onyango 2002",]) +
  geom_line(data=sims,aes(x=DPIDissorTrans,y=meanDiss,group=Moz,col=Moz),linetype=3,alpha=0.7) +
  geom_errorbar(aes(x=DPIDissorTrans,ymin=lowerDiss,ymax=upperDiss),alpha=0.5) +
   geom_point(aes(x=DPIDissorTrans,y=meanDiss,fill=factor(Moz)),shape=21) +
  scale_colour_manual(values=c("royalblue4","dodgerblue")) +
  scale_fill_manual(values=c("royalblue4","dodgerblue")) +
  labs(fill="",title="B") +
   xlab("Days post blood meal") +
  ylab("Proportion of mosquitoes with a disseminated infection") +
 #facet_wrap(~Moz,labeller=label_wrap_gen(width=8))+
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
        ,legend.position=c(0.8,0.3)
  )  


pdf(file="fig_midgutHemocoelModelFits.pdf",width=6,height=4)
grid.arrange(infFitsPlot,dissFitsPlot,ncol=2)
dev.off()