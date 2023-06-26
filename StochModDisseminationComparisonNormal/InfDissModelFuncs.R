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
                       , tl.params=list(epsilon=0.01)) 
  return(data.frame(out))
}


dissSummaryFunc <- function(modelOutput){
  modelOutput$days <- round(modelOutput$time/24,1)  # create a column of .1 days
  # for each simulation (mosquito) 
  lastRunPerDay <- lapply(unique(modelOutput$run),function(y){
    runD <- modelOutput[modelOutput$run %in% y,]  # select single run
    # for each '10th of a day' in each model run select the row that has the last time point and return the results
    maxTimePerDay <- lapply(unique(runD$days),function(z){
      mtemp <- runD[runD$days %in% z,]
      maxT <- mtemp[mtemp$t %in% max(mtemp$t),]
      return(maxT[length(maxT[,1]),])
    })
    maxTimePerDay <- do.call(rbind.data.frame,maxTimePerDay) # return the results for each daily last time point
    return(maxTimePerDay)
  })
  lastRunPerDay <- do.call(rbind.data.frame,lastRunPerDay)
  # for each 10th of a day, see if any rows have Hc > 0 
  tempDiss <- lapply(unique(lastRunPerDay$days),function(a){
    subDat <- lastRunPerDay[lastRunPerDay$days %in% a,]
    HcInf <- length(subDat$Hc[subDat$Hc>0])
    totalSize <- length(subDat$Hc)
    #propDiss <- HcInf/ length(unique(subDat$run))
    return(c(a,HcInf,totalSize))
  })
  tempDiss2 <- do.call(rbind.data.frame,tempDiss)
  names(tempDiss2) <- c("time","numberRunsDisseminated","totalSize")

  return(tempDiss2)
}
