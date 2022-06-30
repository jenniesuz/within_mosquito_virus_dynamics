
library(GillespieSSA)
library(ggplot2)
library(gridExtra)
library(plyr)
library(parallel)
require(sensitivity)

#**************************************************************************

dissSummaryFunc <- function(dissModelOutput){
  dissModelOutput$days <- round(dissModelOutput$t/24,1)  # create a column of .1 days
  
  lastRunPerDay <- lapply(unique(dissModelOutput$run),function(y){
    runD <- dissModelOutput[dissModelOutput$run %in% y,]  # select single run
    maxTimePerDay <- lapply(unique(runD$days),function(z){
      mtemp <- runD[runD$days %in% z,]
      maxT <- mtemp[mtemp$t %in% max(mtemp$t),]
      return(maxT[1,])
    })
    maxTimePerDay <- do.call(rbind.data.frame,maxTimePerDay)
    return(maxTimePerDay)
  })
  lastRunPerDay <- do.call(rbind.data.frame,lastRunPerDay)
  
  tempDiss <- lapply(unique(lastRunPerDay$days),function(a){
    subDat <- lastRunPerDay[lastRunPerDay$days %in% a,]
    HciInf <- length(subDat$Hci[subDat$Hci>0])
    propDiss <- HciInf/ length(unique(subDat$run))
    return(c(a,propDiss))
  })
  tempDiss2 <- do.call(rbind.data.frame,tempDiss)
  names(tempDiss2) <- c("time","proportionDisseminated")
  return(tempDiss2)
}


#************************Varying input virus concentration**********
virusConcs <- c(10^5,10^6,10^7,10^8)

sim.func <- function(x){
  initial <- c(Gv=x*0.003,Mci=0,Mv=0,Hv=0,Hci=0) # initial state account for average size of moz bm
  
  finalTime <- 336
  
  nsims <- 30 #100
  
  parameters <- c(muV = 0.1
                  ,infRate = 10^-7
                  ,prodRate =   10  
                  ,cellSpread = 10^-6
                  ,escapeRate = 0.005  
                  ,cMax = 400   
                  ,hMax = 900             
  )
  
  a <- c("muV*Gv"
         ,"Gv*infRate*(cMax-Mci)"
         ,"cellSpread*Mci*(cMax-Mci)"
         ,"prodRate*Mci"
         ,"muV*Mv"
         ,"escapeRate*Mv"
         ,"prodRate*Hci"
         ,"muV*Hv"
         ,"Hv*infRate*(hMax-Hci)"
  ) # character vector of propensity functions
  
  nu <- cbind(c(-1,0,0,0,0)
              ,c(-1,+1,0,0,0)
              ,c(0,+1,0,0,0)
              ,c(0,0,+1,0,0)
              ,c(0,0,-1,0,0)
              ,c(0,0,-1,+1,0)
              ,c(0,0,0,+1,0)
              ,c(0,0,0,-1,0)
              ,c(0,0,0,0,+1))
  
  repSim <- lapply(1:nsims,function(y){
    
  out <- ssa(initial,a,nu,parameters,tf=finalTime,method="ETL")
  out <- data.frame(out$dat)
     
  maxTime <-  max(out$t)
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
     
    out$run <- y
    
    out$inf <- 0
    if(out$Mv[length(out$Mv)]>0){out$inf <- 1}
    return(out)
  })
  
  repSims <- do.call(rbind.data.frame,repSim)
  repSims$conc <- x
  repSims$MciProp <- repSims$Mci/parameters[6]
  return(repSims)

}

set.seed(123)

doseSim <- mclapply(virusConcs, sim.func)

doseSim2 <- do.call(rbind.data.frame,doseSim)
doseSim2$virusConc <- log10(doseSim2$conc)

#***********Change in numbers of virions in the midgut epi over time****
test <- doseSim2

# calculate time to midgut infection

timeToMG <- lapply(unique(test$conc),function(x){
  temp <- test[test$conc %in% x,]
  minTRuns <- lapply(unique(test$run),function(y){
    temp2 <- temp[temp$run %in% y,]
    temp2 <- temp2[!temp2$Mci %in% 0,]
    if(length(temp2$Mci)>0){
    minT <- temp2$t[min(temp2$Mci)]
    }else{
      minT<- NA
    }
    return(cbind(minT,run=y))
  })
  minTRuns <- do.call(rbind.data.frame,minTRuns)
  minTRuns$conc <- x
  return(minTRuns)
})

timeToMG <- do.call(rbind.data.frame,timeToMG)
timeToMG

ddply(timeToMG,.(conc),summarise,minI=min(minT,na.rm=T),maxI=max(minT,na.rm=T))

#********range in number of virions in the midgut across simulations
testSub <- test[test$conc %in% 1e+07,]

# 24hrs
testSub2 <- testSub[(testSub$t >24 )& (testSub$t<=25),]
byRun <- ddply(testSub2,.(run),summarise,max=max(Mv,na.rm=T))
range(byRun$max)

# 72hrs
testSub2 <- testSub[(testSub$t >72 )& (testSub$t<=73),]
byRun <- ddply(testSub2,.(run),summarise,max=max(Mv,na.rm=T))
range(byRun$max)



#***********Plot change in numbers of virions in the midgut epi over time****

mVPlot <- ggplot(test) +
  geom_line(aes(x=t,y=Mv,group=run),col="grey",lwd=0.25) +
  xlim(0,48) +
  labs(title="A") +
  xlab("Time (hours)") +
  ylab("Number of virions in the midgut epithelium (Mv)") +
  facet_wrap(~virusConc,scales="free",labeller = label_bquote("Number of input virions (Gv):"~.(10^virusConc*0.003))) +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,strip.background = element_rect(colour="white", fill="white")
        ,strip.text = element_text(size = 6)
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=6)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=6)
        ,legend.position ="none"
        # ,legend.title = element_blank()
  )

hVPlot <- ggplot(test) +
  geom_line(aes(x=t,y=Hv,group=run),col="grey",lwd=0.25) +
  xlab("Time (hours)") +
  labs(title="B") +
  ylab("Number of virions in the hemocoel (Hv)") +
  facet_wrap(~virusConc,scales="free",labeller = label_bquote("Input number of virions (Gv):"~.(10^virusConc*0.003))) +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,strip.background = element_rect(colour="white", fill="white")
        ,strip.text = element_text(size = 6)
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=6)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=6)
        ,legend.position ="none"
        # ,legend.title = element_blank()
  )


pdf(file="fig_stochVarSimsConc.pdf",width=4,height=8)
grid.arrange(mVPlot,hVPlot)
dev.off()

#**************************************


#**********first remove simulations where an infection didn't occur**********************
inf <- ddply(doseSim2,.(run,conc),summarise,occurred=sum(inf))  # establish if infection occurred
inf$occurred[inf$occurred>0] <- 1
inf$runConc <- paste(inf$run,inf$conc)
noInf <- inf$runConc[inf$occurred == 0]                         # note run and concs where infection didn't occur

doseSim2$runConc <- paste(doseSim2$run,doseSim2$conc)           # remove these from the dissemination dataset
doseSim2 <- doseSim2[!doseSim2$runConc %in% noInf,]
#*******************************************************************************

diss <- lapply(unique(doseSim2$virusConc),function(x){
  temp <- doseSim2[doseSim2$virusConc %in% x,]   # for each input virus concentration
  tempDiss <-dissSummaryFunc(temp)
  tempDiss$conc <- x
  return(tempDiss)
})

diss2 <- do.call(rbind.data.frame,diss)




#*****************Plot*************
propDiss <- ggplot(diss2) +
  geom_line(aes(x=time,y=proportionDisseminated,col=conc,group=conc)) +
  xlab("Time (days)") +
  ylab("Proportion of simulations with disseminated \n infection given infection in midgut") +
  ylim(0,1) +
  labs(col=expression("Input virus particles (log"[10]*")")) +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=8)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,legend.position =c(0.8,0.5)
  ) 

propDiss

pdf(file="fig_hemocoelVirusConc.pdf",width=5,height=5)
propDiss
dev.off()


#**************************************************************************


















#********************run for limited values of the virus production rate**********************************
gamma <- c(1,10,100,1000)
run <- c(1:length(gamma))
randSnd <- cbind.data.frame(gamma,run)
#*******************************Run simulations in parallel****************************
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
#*************************************************************************************
sensGamma <- parRapply(cl,randSnd,function(y){

  library(GillespieSSA)
  library(plyr)
  library(parallel)
  
  params= c(muV = 0.1
            ,infRate = 10^-7
            ,prodRate =   as.numeric(y[1]) 
            ,cellSpread = 10^-6 
            ,escapeRate = 0.005
            ,cMax = 400   
            ,hMax = 900             
  )

  source("4_sensitivityAnalysisHemocoelSimFunc.R")

  sim <- sim.func(x=10^8,parms=params)
  sim$run <- as.numeric(y[2])
  return(sim)

})

stopCluster(cl)

sensGamma <- do.call(rbind.data.frame,sensGamma)


diss <- lapply(unique(sensGamma$run),function(x){
  temp <- sensGamma[sensGamma$run %in% x,]
  temp$t <- round(temp$t,0)
  tempDiss <- lapply(unique(temp$t),function(y){ # for all runs sample a time point
    subDat <- temp[temp$t %in% y,]               # at that time point total number of runs that
    HciInf <- length(subDat$Hci[subDat$Hci>0])   # have Hci>0
    MciInf <- length(subDat$Mci[subDat$Mci>0])
    if(MciInf==0){propDiss<-0}else{
      propDiss <- HciInf/ MciInf 
    }# calculate proportion 
    return(c(y,propDiss))
  })
  tempDiss2 <- do.call(rbind.data.frame,tempDiss)
  names(tempDiss2) <- c("time","proportionDisseminated")
  tempDiss2$run <- x
  return(tempDiss2)
})

dissGamma <- do.call(rbind.data.frame,diss)
dissGamma <- dissGamma[!dissGamma$time %in% max(dissGamma$time),]
dissGamma2 <- merge(dissGamma,randSnd,by="run",all.x=T)
#******************************Plot these results***********************************
gammaPlot <- ggplot(dissGamma2) +
  geom_line(aes(x=time/24,y=proportionDisseminated,group=as.factor(run),col=gamma)) +
  xlab(expression("Time (days)")) +
  ylab("Proportion of simulations with disseminated infection") +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=10)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=10)
     #  ,legend.position ="none"
     #  ,legend.title = element_blank()
  )

gammaPlot







#****************************Run for limited values of virus decay rate*******************************
decay <- c(1/12,1/9,1/6,1/3)
run <- c(1:length(decay))
randSnd <- cbind.data.frame(decay,run)
#*******************************Run simulations in parallel****************************
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
#*************************************************************************************
sensDecay <- parRapply(cl,randSnd,function(y){
  
  library(GillespieSSA)
  library(plyr)
  library(parallel)
  
  params= c(muV = as.numeric(y[1])
            ,infRate = 10^-7
            ,prodRate =   10 
            ,cellSpread = 10^-6 
            ,escapeRate = 0.005 
            ,cMax = 400   
            ,hMax = 900             
  )
  
  source("4_sensitivityAnalysisHemocoelSimFunc.R")
  
  sim <- sim.func(x=10^8,parms=params)
  sim$run <- y[2]
  return(cbind.data.frame(sim))
  
})

stopCluster(cl)

sensDecay <- do.call(rbind.data.frame,sensDecay)


diss <- lapply(unique(sensDecay$run),function(x){
  temp <- sensDecay[sensDecay$run %in% x,]
  temp$t <- round(temp$t,0)
  tempDiss <- lapply(unique(temp$t),function(y){ # for all runs sample a time point
    subDat <- temp[temp$t %in% y,]               # at that time point total number of runs that
    HciInf <- length(subDat$Hci[subDat$Hci>0])   # have Hci>0
    MciInf <- length(subDat$Mci[subDat$Mci>0])
    if(MciInf==0){propDiss<-0}else{
      propDiss <- HciInf/ MciInf 
    }# calculate proportion 
    return(c(y,propDiss))
  })
  tempDiss2 <- do.call(rbind.data.frame,tempDiss)
  names(tempDiss2) <- c("time","proportionDisseminated")
  tempDiss2$run <- x
  return(tempDiss2)
})
dissDecay <- do.call(rbind.data.frame,diss)
dissDecay <- dissDecay[!dissDecay$time %in% max(dissDecay$time),]
dissDecay2 <- merge(dissDecay,randSnd,by="run",all.x=T)
#******************************Plot these results***********************************
decayPlot <- ggplot(dissDecay2) +
  geom_line(aes(x=time/24,y=proportionDisseminated,group=as.factor(run),col=decay)) +
  xlab(expression("Time (days)")) +
  ylab("Proportion of simulations with disseminated infection") +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=10)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=10)
        #  ,legend.position ="none"
        #  ,legend.title = element_blank()
  )

decayPlot


#*********************************



