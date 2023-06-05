#**************************Function to wrap fitting and simulations for iteration*************
repeatModelRunFunc <- function(x
                               ,dpid=dat$DPIDissorTrans
                               ,conc=dat$Conc.Min
                               ,p1=parms$muV
                               ,p2=parms$infRate1
                               ,p3=parms$prodRate1
                               ,p4=parms$cellSpread1
                               ,p5=parms$escapeRate1
                               ,p6=parms$cMax
                               ,p7=parms$hMax
){
  a<-x
  sim <- sim.func(x=10^conc[1]
                  ,muV=p1
                  ,infRate=p2
                  ,prodRate=p3
                  ,cellSpread=p4
                  ,escapeRate=p5
                  ,cMax=p6
                  ,hMax=p7)
  
  #**********first remove simulations where an infection didn't occur**********************
  inf <- ddply(sim,.(run),summarise,occurred=sum(inf))  # establish if infection occurred
  inf$occurred[inf$occurred>0] <- 1
  noInf <- inf$run[inf$occurred == 0]                         # note run and concs where infection didn't occur
  # remove these from the dissemination dataset
  sim <- sim[!sim$run %in% noInf,]
  #*******************************************************************************
  tempDiss <- dissSummaryFunc(sim)
  
  subSim <- tempDiss[tempDiss$time %in% as.numeric(dpid),]
  subSim$run2 <- x
  return(subSim)
}
#********************************************************************************
