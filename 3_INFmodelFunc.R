library(GillespieSSA)

sim.func <- function(x,parms) with(c(parms), {
  initial <- c(Gv=x*0.003,Mci=0,Mv=0) # initial state # assume input virus concentration then mosquito imbibes c. 3ul
  
  finalTime <- 124
  
  nsims <- 30
  
  a <- c("muV*Gv"
         ,"Gv*infRate*(cMax-Mci)"
         ,"cellSpread*Mci*(cMax-Mci)"
         ,"prodRate*Mci"
         ,"muV*Mv"
         ,"escapeRate*Mv"
  ) # character vector of propensity functions
  #
  nu <- cbind(c(-1,0,0)
              ,c(-1,+1,0)
              ,c(0,+1,0)
              ,c(0,0,+1)
              ,c(0,0,-1)
              ,c(0,0,-1)
  )
  
  repSim <- lapply(1:nsims,function(z){
    out <- ssa(initial,a,nu,parms,tf=finalTime,method="ETL")
    
    dat<-data.frame(out$dat)
    dat$run <- z
    
    dat$inf <- 0
    
    if(dat$Mv[length(dat$Mv)]>0){dat$inf <- 1}
    return(dat)
  })
  
  repSims <- do.call(rbind.data.frame,repSim)
  repSims$conc <- x
  repSims$MciProp <- repSims$Mci/as.numeric(parms[6])
  return(repSims)
 })




