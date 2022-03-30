library(GillespieSSA)

sim.func <- function(x,parms) with(c(parms), {
  initial <- c(Gv=round(x*0.003,0),Mci=0,Mv=0,Hv=0,Hci=0) # initial state # assume input virus concentration then mosquito imbibes c. 3ul

  finalTime <- 336

  nsims <- 30

  a <- c("muV*Gv"
         ,"infRate*Gv*(cMax-Mci)"
         ,"cellSpread*Mci*(cMax-Mci)"
         ,"prodRate*Mci"
         ,"muV*Mv"
         ,"escapeRate*Mv"
         ,"prodRate*Hci"
         ,"muV*Hv"
         ,"Hv*infRate*(hMax-Hci)") # character vector of propensity functions
  #
  nu <- cbind(c(-1,0,0,0,0)
              ,c(-1,+1,0,0,0)
              ,c(0,+1,0,0,0)
              ,c(0,0,+1,0,0)
              ,c(0,0,-1,0,0)
              ,c(0,0,-1,+1,0)
              ,c(0,0,0,+1,0)
              ,c(0,0,0,-1,0)
              ,c(0,0,0,0,+1))

  repSim <- lapply(1:nsims,function(z){
    
   out <- ssa(initial,a,nu,parms,tf=finalTime,method="ETL")
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
        
    out$run <- z

    out$inf <- 0

    if(out$Mv[length(out$Mv)]>0){out$inf <- 1}
    return(out)
  })

  repSims <- do.call(rbind.data.frame,repSim)
  repSims$conc <- x
  return(repSims)
 })




