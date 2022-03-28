library(deSolve)
library(GillespieSSA)
library(ggplot2)
library(gridExtra)
library(grid)
library(plyr)
library(parallel)
require(sensitivity)
library(lhs)

# Sensitivity analysis varying the rate of cell infection and the midgut virus decay rate to simulate the dose-response curves of the 
# probability of midgut infection


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
# saveRDS(sensSims,"sensitivitySubsetBeta140921.rds")
#***************************************************************


#************Run saved sims from above*********
beta <- c(-14,-12,-10,-8,-6,-4)
run <- c(1:length(beta))
randSnd <- cbind.data.frame(beta,run)
sensSims <- readRDS("sensitivitySubsetBeta140921.rds")

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

#pdf(file="fig_midgutPropInfBeta.pdf",width=5,height=5)
#pInfBeta
#dev.off()


#********************************************************************************************
#****************************Run for limited values of virus decay rate*******************************


# #********************run for limited values of beta**********************************
# decay <- c(1/36,1/30,1/24,1/18,1/12,1/6)
# run <- c(1:length(decay))
# randSnd <- cbind.data.frame(decay,run)
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
# saveRDS(sensSims,"sensitivitySubsetDecay140921.rds")




#**************read in saved sims from above***************
decay <- c(1/36,1/30,1/24,1/18,1/12,1/6)
run <- c(1:length(decay))
randSnd <- cbind.data.frame(decay,run)

sensSims <- readRDS("sensitivitySubsetDecay140921.rds")
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

pdf(file="fig_midgutPropInfDecayBeta.pdf",width=5,height=5)
grid.arrange(pInfBeta,pInfDecay,ncol=1
             ,left = textGrob("Proportion of simulations with established midgut infection", rot = 90, vjust = 1,gp=gpar(fontsize=10))
             )
dev.off()











#*****************************Read in from earlier run and plot*********************
sensSims <- readRDS("sensitivity130921.rds")

infDat2 <- read.csv("210913Sensitivity.csv")



#write.csv(infDat2,"210913Sensitivity.csv")
#*********************************************************************

infDat2$denom<-30

mods <- sapply(unique(infDat2$run),function(x){
  temp <- infDat2[infDat2$run %in% x,]
  if(max(temp$propInf)>0){
  plot(temp$virusConc,temp$propInf)
  modTest <- glm(propInf ~ virusConc+beta,family="binomial",data=temp,weights=denom)
  return(cbind(x,exp(coef(modTest)[2])))}else{
    return(cbind(x,NA))
  }
  
})


mods <- cbind.data.frame(mods[1,],mods[2,])
names(mods)<- c("run","oddsInf")
mods$beta <- randSnd$beta
mods$decay <- randSnd$virusDecay

ggplot(mods) +
  geom_point(aes(beta,oddsInf,col=decay)) +
  ylim(0,100)

ggplot(mods) +
  geom_point(aes(decay,oddsInf,col=beta)) +
  ylim(0,100)


runs <- mods$run[mods$oddsInf < 10]

infDat3 <- infDat2[infDat2$run %in% runs,]

ggplot(infDat3) +
  geom_line(aes(x=virusConc,y=propInf,group=as.factor(run),col=as.factor(run))) +
  xlab(expression("Input virus concentration (log"[10]*")")) +
  ylab("Proportion of simulations with established midgut infection") +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=10)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=10)
        ,legend.position ="none"
        # ,legend.title = element_blank()
  )


#****************
# virus decay rate in bloodmeal no impact between average of 5 and 76 hrs
# why?
# what time point does virus infection of midgut cells occur?
# min t where Mci >0

min(sensSims$t[sensSims$Mci>0])

timeTomg <- ddply(sensSims,.(run,conc,sampleNum),summarise,tMci=min(t[which(Mci>0)]))
timeToMg <- timeTomg[!timeTomg$tMci %in% Inf,]
timeToMg <- merge(timeToMg,randSnd,by="run",all.x=T)

range(timeToMg)
max(timeToMg$tMci)

ggplot(timeToMg) +
  geom_point(aes(x=log10(conc),y=tMci/24,group=run,col=beta))

ggplot(timeToMg) +
  geom_point(aes(x=log10(conc),y=tMci/24,group=run,col=virusDecay))

ggplot(timeToMg) +
  geom_point(aes(x=virusDecay,y=tMci,group=conc,col=log10(conc)))
ggplot(timeToMg) +
  geom_point(aes(x=beta,y=tMci,group=conc,col=log10(conc)))
ggplot(timeToMg) +
  geom_point(aes(x=beta,y=tMci,group=conc,col=log10(conc))) +
  ylim(0,24)


ggplot(timeToMg) +
  geom_histogram(aes(tMci)) 


mean(timeToMg$tMci) # mean of 3.3 hrs.

# virus decay
ggplot(sensSims[sensSims$run %in% 1:5,]) +
  geom_point(aes(x=t,y=Gv,group=sampleNum)) +
  facet_wrap(~run)
  





#******************************Generate parameter values****************************

#virusDecay <- c(0.013,0.2)  # from about 3 days to about 5 hours
#beta <- c(-15,-5)

#min <- c(virusDecay[1],beta[1])
#max <- c(virusDecay[2],beta[2])

#params <- c("virusDecay"
#            ,"beta" 
#         )

#params <- cbind.data.frame(params,min,max)

#r <- randomLHS(100,length(params[,1]))
#parmVals <- lapply(1:length(params[,1]),function(x){
#  temp <- params[x,]
#  randomSample <- runif(r[,x],min=temp$min,max=temp$max)

#})
#parmVals <- do.call(cbind.data.frame,parmVals)
#names(parmVals) <- params$params

#randSnd <- parmVals
#randSnd$run <- 1:100

#*******************************Run simulations in parallel****************************

#no_cores <- detectCores() - 1
#cl <- makeCluster(no_cores)

#*************************************************************************************

#sensSims <- parRapply(cl,randSnd,function(y){

#  library(GillespieSSA)
#  library(plyr)
#  library(parallel)
#  virusConcs <- c(10^3,10^3.5,10^4,10^4.5,10^5,10^5.5,10^6,10^6.5,10^7,10^7.5,10^8,10^8.5,10^9,10^9.5,10^10)

#  params= c(muV = as.numeric(y[1])
#            ,infRate = 10^as.numeric(y[2])
#            ,prodRate =   20  
#            ,cellSpread = 10^-4  
#            ,escapeRate = 0.12  
#            ,cMax = 400   
#            ,hMax = 900             
#  )
#saveRDS(params,"params")

#  sim.func <- function(x,parms=params){
#    initial <- c(Gv=x*0.003,Mci=0,Mv=0,Hv=0,Hci=0) # initial state 
#    
#    finalTime <- 124
#    
#    nsims <- 30
#    
#    parameters <- parms#readRDS("params")
#    
#    a <- c("muV*Gv"
#           ,"Gv*infRate*(cMax-Mci)"
#           ,"cellSpread*Mci*(cMax-Mci)"
#           ,"prodRate*Mci"
#           ,"muV*Mv"
#           ,"escapeRate*Mv"
#           ,"prodRate*Hci"
#           ,"muV*Hv"
#           ,"Hv*infRate*(hMax-Hci)"
#    ) # character vector of propensity functions

#    nu <- cbind(c(-1,0,0,0,0)
#                ,c(-1,+1,0,0,0)
#                ,c(0,+1,0,0,0)
#                ,c(0,0,+1,0,0)
#                ,c(0,0,-1,0,0)
#                ,c(0,0,-1,+1,0)
#                ,c(0,0,0,+1,0)
#                ,c(0,0,0,-1,0)
#                ,c(0,0,0,0,+1))

#    repSim <- lapply(1:nsims,function(z){
#      out <- ssa(initial,a,nu,parameters,tf=finalTime,method="ETL")
#      
#      dat<-data.frame(out$dat)
#      dat$run <- z
#      
#      dat$inf <- 0
#      
#      if(dat$Mv[length(dat$Mv)]>0){dat$inf <- 1}
#      return(dat)
#    })
#    
#    repSims <- do.call(rbind.data.frame,repSim)
#    repSims$conc <- x
#    repSims$MciProp <- repSims$Mci/parameters[6]
#    return(repSims)
#  }


#  doseSim <- mclapply(virusConcs, sim.func)
#  doseSim2 <- do.call(rbind.data.frame,doseSim)
#  run <- y[3]
#  return(cbind.data.frame(doseSim2,run))

#})

#stopCluster(cl)

#sensSims <- do.call(rbind.data.frame,sensSims)
#names(sensSims)[7]<- "sampleNum"
#saveRDS(sensSims,"sensitivity130921.rds")

#infDat <- ddply(sensSims,.(run,conc,sampleNum),summarise,sumRuns=sum(inf))
#infDat$sumRuns[infDat$sumRuns>0] <-1
#infDat2 <- ddply(infDat,.(run,conc),summarise,propInf=sum(sumRuns)/max(sampleNum))
#infDat2$virusConc <- log10(infDat2$conc)

#infDat2 <- merge(infDat2,randSnd,by="run",all.x=T)
#*******************************
