source(here::here("StochModInfectionComparisonBeta//INFmodelFunc.R"))
library(parallel)


repeatInfModel <- function(x=virus_params()
    ,startingVirus=10^competenceDat$ConcMax[competenceDat$Moz %in% "Ae. aegypti"]
){

  simEachConc <- tryCatch ({ 
    lapply(startingVirus,function(conc){
      repSim <- lapply(1:30,function(y){
        out <- infectionModel(startingVirus=conc
                              ,muV = as.numeric(x[1])
                              ,infRate = as.numeric(x[2])
                              ,prodRate = as.numeric(x[3])    
                              ,cellSpread = as.numeric(x[4])  
                              ,escapeRate = as.numeric(x[5])
                              ,cMax = as.numeric(x[6]) )
        dat <- data.frame(out)
        dat$run <- y
        dat$inf <- 0
        if(dat$Mv[length(dat$Mv)]>0){dat$inf<-1}
        return(dat[1,c("run","inf")])  # just return the run, and whether infection was established
      })
      
      
      repSims <- do.call(rbind.data.frame,repSim)
      return(data.frame("num"=sum(repSims$inf)
               ,"denom"=length(repSims$inf)
               ,"conc"=log10(conc)))
      
    }) 
    
  } , error = function(e) {
    print(e)
    return(c(NA,NA,NA))
  }
  )
  
  if(is.na(simEachConc[[1]][1])==F){  
    simDat <- do.call(rbind.data.frame,simEachConc)
    names(simDat) <- c("num","denom","conc")
    return(simDat)
  }else{
    return(c(NA,NA,NA))
  }
}

#startTime <- Sys.time()
#test <- repeatInfModel()
#endTime <- Sys.time()
#endTime - startTime