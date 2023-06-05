
library(parallel)

virus_params <- function(   muV = 0.1
                            ,muVh=0.2
                            ,infRate = 10^-8
                            ,prodRate = 10            # note these parameters don't matter for infection
                            ,cellSpread = 10^-4        
                            ,escapeRate = 10^-5         #
                            ,cMax = 400  
                            ,hMax = 400
)
  return(as.list(environment()))


repeatModel <- function(x=virus_params()
    ,startingVirus=10^competenceDat$ConcMax[competenceDat$Moz %in% "Ae. aegypti"]
){

#  simEachConc <- #tryCatch ({ 
    simEachConc <-    lapply(startingVirus,function(conc){
      repSim <- lapply(1:30,function(y){
        out <- infDissModel(startingVirus=conc
                              ,muV = as.numeric(x[1])
                              ,infRate = as.numeric(x[2])
                              ,prodRate = as.numeric(x[3])    
                              ,cellSpread = as.numeric(x[4])  
                              ,escapeRate = as.numeric(x[5])
                              ,cMax = as.numeric(x[6]) 
                              ,hMax = as.numeric(x[7])
        )
        dat <- data.frame(out)
        dat$run <- y
        dat$inf <- 0
        if(dat$Mc[length(dat$Mc)]>0){dat$inf<-1}
        return(dat)
  })
      
      
      repSims <- do.call(rbind.data.frame,repSim)
      return(repSims)
      
    }) 
    
 # } , error = function(e) {
#    print(e)
#    return(c(NA,NA,NA))
#  }
#  )
  
 # if(is.na(simEachConc[[1]][1])==F){  
#    simDat <- do.call(rbind.data.frame,simEachConc)
#    names(simDat) <- c("num","denom","conc")
#    return(simDat)
#  }else{
#    return(c(NA,NA,NA))
 # }
}

startTime <- Sys.time()
test <- repeatModel(startingVirus=10^6)
endTime <- Sys.time()
endTime - startTime
