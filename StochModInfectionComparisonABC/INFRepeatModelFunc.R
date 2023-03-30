source(here::here("StochModInfectionComparisonABC//INFmodelFunc.R"))

repeatInfModel <- function(nsims=30,sv=10^8){

  cl <- makeCluster(detectCores()-1)
  clusterEvalQ(cl, {library(adaptivetau)})
  environment(infectionModel) <- .GlobalEnv
  clusterExport(cl, varlist=c("infectionModel","sv","nsims"),
                envir=environment())
  
  repSim <- parLapply(cl,1:nsims,function(y){
    out <- infectionModel(startingVirus=sv)
     dat<-data.frame(out)
    dat$run <- y
    dat$inf <- 0
    if(dat$Mv[length(dat$Mv)]>0){dat$inf<-1}
    return(dat)
  })

  stopCluster()
  
  repSims <- do.call(rbind.data.frame,repSim)

  repSims$conc <- x
  repSims$MciProp <- repSims$Mci/as.numeric(params[6])
  return(repSims)
}

