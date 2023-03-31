library(deSolve)
library(GillespieSSA)
library(ggplot2)
library(gridExtra)
library(grid)
library(parallel)
library(plyr)

# Sensitivity analysis varying the rate of cell infection and the midgut virus decay rate to simulate the dose-response curves of the 
# probability of midgut infection

# the below code has been run and outputs saved as rds files.

#********************run for limited values of beta**********************************
# beta <- c(-14,-12,-10,-8,-6,-4)
# run <- c(1:length(beta))
# randSnd <- cbind.data.frame(beta,run)
# #*******************************Run simulations in parallel****************************
# no_cores <- detectCores() - 1
# cl <- makeCluster(no_cores)
# #*************************************************************************************
# sensSims2 <- parRapply(cl,randSnd,function(y){
# 
#   library(GillespieSSA)
#   library(plyr)
#   library(parallel)
#   virusConcs <- c(10^3,10^3.5,10^4,10^4.5,10^5,10^5.5,10^6,10^6.5,10^7,10^7.5,10^8,10^8.5,10^9,10^9.5,10^10)
# 
#   params= c(muV = 0.02 # 2 days
#             ,infRate = 10^as.numeric(y[1])
#             ,prodRate =   20
#             ,cellSpread = 10^-4
#             ,escapeRate = 0.12
#             ,cMax = 400
#             ,hMax = 900
#   )
# #saveRDS(params,"params")
# 
#   sim.func <- function(x,parms=params){
#     initial <- c(Gv=x*0.003,Mci=0,Mv=0,Hv=0,Hci=0) # initial state
# 
#     finalTime <- 124
# 
#     nsims <- 30
# 
#     parameters <- parms#readRDS("params")
# 
#     a <- c("muV*Gv"
#            ,"Gv*infRate*(cMax-Mci)"
#            ,"cellSpread*Mci*(cMax-Mci)"
#            ,"prodRate*Mci"
#            ,"muV*Mv"
#            ,"escapeRate*Mv"
#            ,"prodRate*Hci"
#            ,"muV*Hv"
#            ,"Hv*infRate*(hMax-Hci)"
#     ) # character vector of propensity functions
# 
#     nu <- cbind(c(-1,0,0,0,0)
#                 ,c(-1,+1,0,0,0)
#                 ,c(0,+1,0,0,0)
#                 ,c(0,0,+1,0,0)
#                 ,c(0,0,-1,0,0)
#                 ,c(0,0,-1,+1,0)
#                 ,c(0,0,0,+1,0)
#                 ,c(0,0,0,-1,0)
#                 ,c(0,0,0,0,+1))
# 
#     repSim <- lapply(1:nsims,function(z){
#       out <- ssa(initial,a,nu,parameters,tf=finalTime,method="ETL")
# 
#       dat<-data.frame(out$dat)
#       dat$run <- z
# 
#       dat$inf <- 0
# 
#       if(dat$Mv[length(dat$Mv)]>0){dat$inf <- 1}
#       return(dat)
#     })
# 
#     repSims <- do.call(rbind.data.frame,repSim)
#     repSims$conc <- x
#     repSims$MciProp <- repSims$Mci/parameters[6]
#     return(repSims)
#   }
# 
# 
#   doseSim <- mclapply(virusConcs, sim.func)
#   doseSim2 <- do.call(rbind.data.frame,doseSim)
#   run <- y[2]
#   return(cbind.data.frame(doseSim2,run))
# 
# })
# 
# stopCluster(cl)
# 
# sensSims <- do.call(rbind.data.frame,sensSims2)
# names(sensSims)[7]<- "sampleNum"
# saveRDS(sensSims,"sensitivityInfBeta.rds")
#***************************************************************

#************Run saved sims from above*********
beta <- c(-14,-12,-10,-8,-6,-4)
run <- c(1:length(beta))
randSnd <- cbind.data.frame(beta,run)
sensSims <- readRDS(here::here("rdsOutputs//sensitivityInfBeta.rds"))

infDat <- ddply(sensSims,.(run,conc,sampleNum),summarise,sumRuns=sum(inf))
infDat$sumRuns[infDat$sumRuns>0] <-1
infDat2 <- ddply(infDat,.(run,conc),summarise,propInf=sum(sumRuns)/max(sampleNum))
infDat2$virusConc <- log10(infDat2$conc)

infDat2 <- merge(infDat2,randSnd,by="run",all.x=T)
infDat2$denom <- 30
#******************************Plot these results***********************************


mods <- lapply(beta,function(x){
  temp <- infDat2[infDat2$beta %in% x,]
  modTest <- glm(propInf ~ virusConc,family="binomial",data=temp,weights=denom)
  newDat <- data.frame("virusConc"=seq(2,10,0.1))
  preds <- predict(modTest,newdata=newDat,type="response")
  return(cbind.data.frame(preds,x,newDat))
})

mods <- do.call(rbind.data.frame,mods)

pInfBeta <- ggplot(infDat2) +
  geom_point(aes(x=virusConc,y=propInf,group=as.factor(run),col=beta)) +
  geom_line(data=mods,aes(x=virusConc,y=preds,group=as.factor(x),col=x)) +
  xlim(4,10) +
  labs(col=expression(paste(italic(beta)," (log"[10]*")"))
       ,title="A") +
  xlab(" ") +
  ylab(" ") +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=10)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=10)
        #,legend.position =c(0.9,0.7)
        # ,legend.title = element_blank()
  )

pInfBeta

#********************************************************************************************
#****************************Run for limited values of virus decay rate*******************************


# #********************run for limited values of beta**********************************
# decay <- c(1/36,1/30,1/24,1/18,1/12,1/6)
# run <- c(1:length(decay))
# randSnd <- cbind.data.frame(decay,run)
#*******************************Run simulations in parallel****************************
# no_cores <- detectCores() - 1
# cl <- makeCluster(no_cores)
#*************************************************************************************
# sensSims2 <- parRapply(cl,randSnd,function(y){
# 
#   library(GillespieSSA)
#   library(plyr)
#   library(parallel)
#   virusConcs <- c(10^3,10^3.5,10^4,10^4.5,10^5,10^5.5,10^6,10^6.5,10^7,10^7.5,10^8,10^8.5,10^9,10^9.5,10^10)
# 
#   params= c(muV = as.numeric(y[1])
#             ,infRate = 10^-8
#             ,prodRate =   20
#             ,cellSpread = 10^-4
#             ,escapeRate = 0.12
#             ,cMax = 400
#             ,hMax = 900
#   )
#   #saveRDS(params,"params")
# 
#   sim.func <- function(x,parms=params){
#     initial <- c(Gv=x*0.003,Mci=0,Mv=0,Hv=0,Hci=0) # initial state
# 
#     finalTime <- 124
# 
#     nsims <- 30
# 
#     parameters <- parms#readRDS("params")
# 
#     a <- c("muV*Gv"
#            ,"Gv*infRate*(cMax-Mci)"
#            ,"cellSpread*Mci*(cMax-Mci)"
#            ,"prodRate*Mci"
#            ,"muV*Mv"
#            ,"escapeRate*Mv"
#            ,"prodRate*Hci"
#            ,"muV*Hv"
#            ,"Hv*infRate*(hMax-Hci)"
#     ) # character vector of propensity functions
# 
#     nu <- cbind(c(-1,0,0,0,0)
#                 ,c(-1,+1,0,0,0)
#                 ,c(0,+1,0,0,0)
#                 ,c(0,0,+1,0,0)
#                 ,c(0,0,-1,0,0)
#                 ,c(0,0,-1,+1,0)
#                 ,c(0,0,0,+1,0)
#                 ,c(0,0,0,-1,0)
#                 ,c(0,0,0,0,+1))
# 
#     repSim <- lapply(1:nsims,function(z){
#       out <- ssa(initial,a,nu,parameters,tf=finalTime,method="ETL")
# 
#       dat<-data.frame(out$dat)
#       dat$run <- z
# 
#       dat$inf <- 0
# 
#       if(dat$Mv[length(dat$Mv)]>0){dat$inf <- 1}
#       return(dat)
#     })
# 
#     repSims <- do.call(rbind.data.frame,repSim)
#     repSims$conc <- x
#     repSims$MciProp <- repSims$Mci/parameters[6]
#     return(repSims)
#   }
# 
# 
#   doseSim <- mclapply(virusConcs, sim.func)
#   doseSim2 <- do.call(rbind.data.frame,doseSim)
#   run <- y[2]
#   return(cbind.data.frame(doseSim2,run))
# 
# })
# 
# stopCluster(cl)
# 
# sensSims <- do.call(rbind.data.frame,sensSims2)
# names(sensSims)[7]<- "sampleNum"
# saveRDS(sensSims,"sensitivityInfDecay.rds")




#**************read in saved sims from above***************
decay <- c(1/36,1/30,1/24,1/18,1/12,1/6)
run <- c(1:length(decay))
randSnd <- cbind.data.frame(decay,run)

sensSims <- readRDS(here::here("rdsOutputs//sensitivityInfDecay.rds"))
infDat <- ddply(sensSims,.(run,conc,sampleNum),summarise,sumRuns=sum(inf))
infDat$sumRuns[infDat$sumRuns>0] <-1
infDat2 <- ddply(infDat,.(run,conc),summarise,propInf=sum(sumRuns)/max(sampleNum))
infDat2$virusConc <- log10(infDat2$conc)

infDat2 <- merge(infDat2,randSnd,by="run",all.x=T)
infDat2$denom <- 30
#******************************Plot these results***********************************

mods <- lapply(decay,function(x){
  temp <- infDat2[infDat2$decay %in% x,]
  modTest <- glm(propInf ~ virusConc,family="binomial",data=temp,weights=denom)
  newDat <- data.frame("virusConc"=seq(2,10,0.1))
  preds <- predict(modTest,newdata=newDat,type="response")
  return(cbind.data.frame(preds,x,newDat))
})

mods <- do.call(rbind.data.frame,mods)

pInfDecay <- ggplot(infDat2) +
  geom_point(aes(x=virusConc,y=propInf,group=as.factor(run),col=decay)) +
  geom_line(data=mods,aes(x=virusConc,y=preds,group=as.factor(x),col=x)) +
  xlim(4,10) +
  labs(col=expression(italic(mu["V"]))
       ,title="B") +
  xlab(expression("Input virus concentration (log"[10]*")  (Gv)")) +
  ylab(" ") +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=10)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=10)
    #    ,legend.position =c(9,0.5)
     #  ,legend.title = element_blank()
  )

pInfDecay

#pdf(file="fig_midgutPropInfDecayBeta.pdf",width=5,height=5)
grid.arrange(pInfBeta,pInfDecay,ncol=1
             ,left = textGrob("Proportion of simulations with established midgut infection", rot = 90, vjust = 1,gp=gpar(fontsize=10))
             )
#dev.off()





