library(adaptivetau)

infDissModel <- function(startingVirus
                           ,muV = 0.1
                           ,infRate = 10^-7.5
                           ,prodRate =   1    
                           ,cellSpread = 10^-4  
                           ,escapeRate = 0.05
                           ,cMax = 400 
                           ,hMax = 400
){
  
  
  params <- list(muV = muV
                 ,probInf = infRate
                 ,prodRate =   prodRate 
                 ,cellSpread = cellSpread
                 ,escapeRate = escapeRate    
                 ,cMax = cMax
                 ,hMax = hMax
  )
  
  transitions <- list(c(Gv = -1,Mc = +1)
                      ,c(Gv = -1)
                      ,c(Mc = +1)
                      ,c(Mv = +1)
                      ,c(Mv = -1)
                      ,c(Mv = -1, Hv = +1)
                      ,c(Hc = +1)
                      ,c(Hv = +1)
                      ,c(Hv = -1)
  )
  
  
  lvrates <- function(y,params,t){
    return( c(y["Gv"]*params$probInf*(params$cMax-y["Mc"])
              ,params$muV*y["Gv"]
              ,params$cellSpread*y["Mc"]*(params$cMax-y["Mc"])
              ,params$prodRate*y["Mc"]
              ,params$muV*y["Mv"]
              ,params$escapeRate*y["Mv"]
              ,y["Hv"]*params$probInf*(params$hMax-y["Hc"])
              ,params$prodRate*y["Hc"]
              ,params$muV*y["Hv"]
    )
    
    )
  }
  
  out <- ssa.adaptivetau(c(Gv = round(startingVirus*0.003,0), Mc = 0, Mv = 0, Hc = 0, Hv = 0),
                       transitions, lvrates, params, tf=400
                       , tl.params=list(epsilon=0.0005)) 
  return(data.frame(out))
}


dissSummaryFunc <- function(modelOutput){
  modelOutput$days <- round(modelOutput$time/24,1)  # create a column of .1 days
  
  lastRunPerDay <- lapply(unique(modelOutput$run),function(y){
    runD <- modelOutput[modelOutput$run %in% y,]  # select single run
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
    HvInf <- length(subDat$Hc[subDat$Hc>0])
    propDiss <- HvInf/ length(unique(subDat$run))
    return(c(a,propDiss))
  })
  tempDiss2 <- do.call(rbind.data.frame,tempDiss)
  names(tempDiss2) <- c("time","proportionDisseminated")

  return(tempDiss2)
}

