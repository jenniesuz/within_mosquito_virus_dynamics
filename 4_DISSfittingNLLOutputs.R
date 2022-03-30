library(gridExtra)
library(ggplot2)
library(plyr)
library(ggfan)
library(binom)

# read in model fits fits
albF <- readRDS("albDissFits220322.rds")
aegF <- readRDS("aegDissFits220322.rds")


# function to change to dataframe and then bind together
listToDat <- function(output=alb, species="alb"){
 params <- lapply(1:length(output),function(x){
  temp <- output[[x]]
  if(is.list(temp)==T){
  temp <- temp[[1]][1]
  temp <- do.call(c,temp)
  names(temp) <- c("prodRate","cellSpread","escapeRate")
  return(exp(temp))
  }else{
    temp <- c(NA,NA,NA)
    names(temp) <- c("prodRate","cellSpread","escapeRate")
    return(temp)
  }
})
 params <- do.call(rbind.data.frame,params)
 names(params) <- c("prodRate","cellSpread","escapeRate")
 params$species=species
return(params)
}

alb <- listToDat(albF)
aeg <- listToDat(output=aegF,species="aeg")

params <- rbind.data.frame(alb,aeg)

params$species[params$species %in% "aeg"] <- "Ae. aegypti"
params$species[params$species %in% "alb"] <- "Ae. albopictus"


pRPlot <- ggplot(params,aes(x=prodRate,fill=species)) + 
  geom_histogram(position="identity",alpha=0.3) +
  scale_fill_manual(values=c("royalblue4","dodgerblue")) +
  scale_x_log10() +
  xlab("Production rate") +
  ylab("Number of fits") +
  labs(fill="") +
  ylab("Number of fits") +
  theme_bw() +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=10)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=10)
        ,legend.position="none"
  )

csPlot <- ggplot(params,aes(x=cellSpread,fill=species)) + 
  geom_histogram(position="identity", alpha = 0.3) +
  scale_fill_manual(values=c("royalblue4","dodgerblue")) +
  scale_x_log10() +
  xlab("Cell spread") +
  ylab("Number of fits") +
  labs(fill="") +
  ylab("Number of fits") +
  theme_bw() +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=10)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=10)
        ,legend.position="none"
  )

erPlot <- ggplot(params,aes(x=escapeRate,fill=species)) + 
  geom_histogram(position="identity", alpha = 0.3) +
  scale_fill_manual(values=c("royalblue4","dodgerblue")) +
  scale_x_log10() +
  xlab("Escape rate") +
  ylab("Number of fits") +
  labs(fill="") +
  ylab("Number of fits") +
  theme_bw() +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=10)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=10)
  )

pdf(file="fig_midgutParameterEstimates.pdf",width=5,height=3)
grid.arrange(pRPlot,csPlot,erPlot,ncol=2)
dev.off()


ddply(params,.(species),summarise,meanCR=mean(contactRate),sdCR=sd(contactRate))
ddply(params,.(species),summarise,meanCR=mean(prodRate),sdCR=sd(prodRate))
ddply(params,.(species),summarise,meanCR=mean(cellSpread),sdCR=sd(cellSpread))
ddply(params,.(species),summarise,meanCR=mean(escapeRate),sdCR=sd(escapeRate))









#***************Simulations**********************
aeg <- readRDS("aegSims.rds")
names(aeg)[7] <- "subRun"

sumAeg <- lapply(unique(aeg$run),function(y){
temp <- aeg[aeg$run %in% y,]  

  tempDiss <- lapply(unique(temp$t),function(z){
   subDat <- temp[temp$t %in% z,]
   MciInf <- length(subDat$Mci[subDat$Mci>0])
   HciInf <- length(subDat$Hci[subDat$Hci>0])
   propDiss <- HciInf/ MciInf
  return(c(z,propDiss))
  }) 
  tempDiss1 <- do.call(rbind.data.frame,tempDiss)
  names(tempDiss1) <- c("time","proportionDisseminated")
  tempDiss1 <- tempDiss1[!tempDiss1$t %in% 0,]
  tempDiss1$days <- round(tempDiss1$time/24,4)
  tempDiss1 <- tempDiss1[,]
  tempDiss1$run <- y
return(tempDiss1)
})

sumAeg2 <- do.call(rbind.data.frame,sumAeg)

ggplot(sumAeg2) +
  geom_line(aes(x=days,y=proportionDisseminated,group=run))

write.csv(sumAeg2,"aegDissSummary220323.csv")


sumAeg <- read.csv("aegDissSummary220323.csv")
names(sumAeg2)[3] <- "DPIDissorTrans"
names(sumAeg2)[2] <- "meanDiss"

ggplot(onyango) +
  geom_line(data=sumAeg2,aes(x=DPIDissorTrans,y=meanDiss,group=run),linetype=3,alpha=0.7,col="royalblue4") +
  geom_errorbar(aes(x=DPIDissorTrans,ymin=lowerDiss,ymax=upperDiss),alpha=0.5) +
  geom_point(aes(x=DPIDissorTrans,y=meanDiss,fill=factor(Moz)),shape=21) +
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













#**********8

alb <- read.csv("albDissSimsSummary.csv")
aeg <- read.csv("aegDissSimsSummary.csv")

alb$species <- ("Ae. albopictus")
aeg$species <- ("Ae. aegypti")


sims <- rbind.data.frame(alb,aeg) 

#******************data to fit to********************
dat <- read.csv("datCiotaOnyango.csv")

bin1 <- binom.confint(x=dat$NumInf,n=dat$ITotal,method="exact")
bin2 <- binom.confint(x=dat$NumDiss,n=dat$DTotal,method="exact")

dat$meanInf <- bin1$mean
dat$lowerInf <- bin1$lower
dat$upperInf <- bin1$upper

dat$meanDiss <- bin2$mean
dat$lowerDiss <- bin2$lower
dat$upperDiss <- bin2$upper


onyango <- dat[dat$Ref %in% "Onyango 2020",]

#****************************************************
names(sims)[2] <- "DPIDissorTrans"
names(sims)[3] <- "meanDiss"
names(sims)[5] <- "Moz"
sims$DPIDissorTrans <- sims$DPIDissorTrans/24

#pdf(file="fig_hemocoelModelFits.pdf",width=5,height=3)
dissFitsPlot <- ggplot(onyango) +
  geom_line(data=sims[sims$Moz %in% "Ae. aegypti",],aes(x=DPIDissorTrans,y=meanDiss,group=run),linetype=3,alpha=0.7,col="royalblue4") +
  geom_line(data=sims[sims$Moz %in% "Ae. albopictus",],aes(x=DPIDissorTrans,y=meanDiss,group=run),linetype=3,alpha=0.7,col="dodgerblue") +
  geom_errorbar(aes(x=DPIDissorTrans,ymin=lowerDiss,ymax=upperDiss),alpha=0.5) +
   geom_point(aes(x=DPIDissorTrans,y=meanDiss,fill=factor(Moz)),shape=21) +
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
#dev.off()
pdf(file="fig_midgutHemocoelModelFits.pdf",width=6,height=4)
grid.arrange(infFitsPlot,dissFitsPlot,ncol=2)
dev.off()