library(MuMIn)
library(here)
library(binom)
library(ggplot2)
library(parallel)
library(adaptivetau)
library(mvtnorm)

# Infection comparison of model fits with observed data
#****model code****
source(here("StochModInfectionComparisonNormal//INFmodelFunc.R"))

source(here("StochModInfectionComparisonNormal//INFRepeatModelFunc.R"))

#******************data to fit to********************
competenceDat <- read.csv(here("StochModInfectionComparisonNormal//datCiotaOnyango.csv"))
competenceDat <- competenceDat[competenceDat$Ref %in% "Ciota 2017",]
bin <- binom.confint(x=competenceDat$NumInf,n=competenceDat$ITotal,method="exact")
competenceDat$meanInf <- bin$mean
competenceDat$lowerInf <- bin$lower
competenceDat$upperInf <- bin$upper
#****************************************************

virus_params <- function(   muV = 0.09676457 
                            ,infRate1 = 1.284587e-09 
                            ,prodRate1 = 1            # note these parameters don't matter for infection
                            ,cellSpread1 = 10^-4        
                            ,escapeRate1 = 0.05         #
                            ,cMax = 400                 #
)
  return(as.list(environment()))


# generate curves for aegypti and albopictus
# using fitted model parameters
#***********aegypti*******
parmsAeg <- virus_params()
nSims <- 30

cl <- makeCluster(detectCores()-1)
clusterEvalQ(cl, {library(adaptivetau)})
environment(repeatInfModel) <- .GlobalEnv
clusterExport(cl, varlist=c("nSims","infectionModel","repeatInfModel","parmsAeg"),
              envir=environment())

simsAeg <- parLapply(cl,1:nSims,function(y){
  #set.seed(y)
  s <- repeatInfModel(x=parmsAeg,startingVirus=10^seq(4,10,length=50))
  s$run <- y
  names(s) <-c("num","denom","conc","run")
  return(s)
})
stopCluster(cl)

simsAeg <- do.call(rbind,simsAeg) 

#**********albopictus*******
virus_params <- function(   muV = 0.09676457 
                            ,infRate1 = 10^-7.519721
                            ,prodRate1 = 1            # note these parameters don't matter for infection
                            ,cellSpread1 = 10^-4        
                            ,escapeRate1 = 0.05         #
                            ,cMax = 400                 #
)
  return(as.list(environment()))
parmsAlb <- virus_params()#infRate=10^-7.519721 )
nSims <- 30

cl <- makeCluster(detectCores()-1)
clusterEvalQ(cl, {library(adaptivetau)})
environment(repeatInfModel) <- .GlobalEnv
clusterExport(cl, varlist=c("nSims","infectionModel","repeatInfModel","parmsAlb"),
              envir=environment())

simsAlb <- parLapply(cl,1:nSims,function(y){
  #set.seed(y)
  s <- repeatInfModel(x=parmsAlb,startingVirus=10^seq(4,8,length=50))
  s$run <- y
  names(s) <-c("num","denom","conc","run")
  return(s)
})
stopCluster(cl)

simsAlb <- do.call(rbind,simsAlb) 
#*********



infPlot <- ggplot(competenceDat) +
  geom_point(aes(x=ConcMax,y=meanInf,col=Moz)) +
  geom_errorbar(aes(x=ConcMax,ymin=lowerInf,ymax=upperInf)) +
  scale_colour_manual(values=c("royalblue4","dodgerblue"),guide="none") +
  geom_line(data=simsAeg,aes(x=conc,y=num/denom,group=run)
            ,linetype = 3
            ,alpha=0.5
            ,col="royalblue4",) +
  geom_line(data=simsAlb,aes(x=conc,y=num/denom,group=run)
            ,linetype = 3
            ,alpha=0.5
            ,col="dodgerblue",) +
  xlab(expression(paste("Virus concentration (log"[10]," PFU/ml)"))) +
  labs(fill="",title="A") +
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

infPlot





#***********************Dissemination****************************
#******************data to fit to********************
competenceDat <- read.csv(here("StochModInfectionComparisonNormal//datCiotaOnyango.csv"))
competenceDat <- competenceDat[competenceDat$Ref %in% "Onyango 2020",]
bin <- binom.confint(x=competenceDat$NumDiss,n=competenceDat$NumInf,method="exact")
competenceDat$meanDiss <- bin$mean
competenceDat$lowerDiss <- bin$lower
competenceDat$upperDiss <- bin$upper

datAeg <- competenceDat[competenceDat$Moz %in% "Ae. aegypti",]
datAlb <- competenceDat[competenceDat$Moz %in% "Ae. albopictus",]
#****************************************************
#****model code****
source(here("StochModDisseminationComparisonNormal//InfDissModelFuncs.R"))
source(here("StochModDisseminationComparisonNormal//InfDissRepeatModelFunc.R"))

virus_params <- function(   muV = 0.09676457 
                            ,infRate = 1.284587e-09 
                            ,prodRate = 100#170         #
                            ,cellSpread = 0.005 # 0.00005       
                            ,escapeRate = 0.005 # 0.25#0.005  #
                            ,cMax = 400  
                            ,hMax = 900             #
)
  return(as.list(environment()))

#*******************ae. aeg**********************

#*
parmsAeg <- virus_params(infRate = 1.284587e-09,prodRate=5.736042e+00,cellSpread=8.939190e-04,escapeRate=3.187694e-01 )
parmsAlb <- virus_params(infRate= 10^-7.519721,prodRate=1.742910e+01,cellSpread= 7.508513e-05,escapeRate=9.880848e-04 )
nSims <- 30
#*
cl <- makeCluster(detectCores()-1)
clusterEvalQ(cl, {library(adaptivetau)})
environment(repeatModel) <- .GlobalEnv
clusterExport(cl, varlist=c("nSims","infDissModel","dissSummaryFunc" ,"repeatModel","datAeg","parmsAeg"),
              envir=environment())

simsAeg <- parLapply(cl,1:nSims,function(y){
  #set.seed(y)
  s <- repeatModel(x=parmsAeg,startingVirus=10^datAeg$ConcMax[1])
  s <- do.call(rbind,s)
  # check for runs that didn't work
  # get rid of runs where midgut infection didn't occur
  s <- s[s$inf %in% 1,]   # make sure calculating the probability of dissemination given infection
  calcDisseminations <- dissSummaryFunc(s)
  calcDisseminations$run <- y
  return(calcDisseminations)
})
stopCluster(cl)

simsAeg <- do.call(rbind,simsAeg) 


#*******************ae. alb**********************
cl <- makeCluster(detectCores()-1)
clusterEvalQ(cl, {library(adaptivetau)})
environment(repeatModel) <- .GlobalEnv
clusterExport(cl, varlist=c("nSims","infDissModel","dissSummaryFunc" ,"repeatModel","datAlb","parmsAlb"),
              envir=environment())

simsAlb <- parLapply(cl,1:nSims,function(y){
  #set.seed(y)
  s <- repeatModel(x=parmsAlb,startingVirus=10^datAlb$ConcMax[1])
  s <- do.call(rbind,s)
  # check for runs that didn't work
  # get rid of runs where midgut infection didn't occur
  s <- s[s$inf %in% 1,]   # make sure calculating the probability of dissemination given infection
  calcDisseminations <- dissSummaryFunc(s)
  calcDisseminations$run <- y
  return(calcDisseminations)
})
stopCluster(cl)

simsAlb <- do.call(rbind,simsAlb) 



dissPlot <- ggplot(competenceDat) +
  geom_line(data=simsAeg,aes(x=time,y=numberRunsDisseminated/totalSize,group=run),linetype=3,alpha=0.7,col="royalblue4") +
  geom_line(data=simsAlb,aes(x=time,y=numberRunsDisseminated/totalSize,group=run),linetype=3,alpha=0.7,col="dodgerblue") +
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

library(grid)
library(gridExtra)

pdf(file="fig_midgutHemocoelmodelfits.pdf",width=5,height=3.5)
grid.arrange(infPlot,dissPlot,ncol=2
             ,left = textGrob("P ", rot = 90, vjust = 1,gp=gpar(fontsize=10))
)
dev.off()
