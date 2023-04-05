source(here::here("StochModInfectionComparisonABC//INFmodelFunc.R"))

repeatInfModel <- function(x=c(s=0.02,par1=0.03,par2=10^-8)
                           ,virusConcs=10^competenceDat$ConcMax[competenceDat$Moz %in% "Ae. aegypti"]
                           ){
  set.seed(x[1])
  simEachConc <- tryCatch ({ 
    lapply(virusConcs,function(conc){

   repSim <- lapply(cl,1:30,function(y){
      out <- infectionModel(startingVirus=conc,par1=x[2],par2=x[3])
       dat<-data.frame(out)
      dat$run <- y
      dat$inf <- 0
      if(dat$Mv[length(dat$Mv)]>0){dat$inf<-1}
      # return(dat)   return the whole timeseries
      return(dat[1,c("run","inf")])  # just return the run,  and whether infection was established
    })

    repSims <- do.call(rbind.data.frame,repSim)
    return(c("prop"=sum(repSims$inf)/(length(repSims$inf))
             ,"num"=sum(repSims$inf)
             ,"denom"=length(repSims$inf)
             ,"conc"=log10(conc)))
  
    }) 
    
    } , error = function(e) {
          print(e)
          return(c(NA,NA,NA,NA))
    }
  )
  
  if(is.na(simEachConc[1])==F){  
  simDat <- do.call(rbind.data.frame,simEachConc)
  names(simDat) <- c("prop","num","denom","conc")
  
  
  modGLM <- glm(matrix(c(num,denom-num),ncol=2) ~ log10(conc)
                ,family="binomial"
                ,data=simDat)
  return(as.numeric(coef(modGLM)))
  }else{
    return(c(-9999,-9999))
  }
}


#test <- repeatInfModel(nsims=30
#                       ,virusConcs=10^competenceDat$ConcMax[competenceDat$Moz %in% "Ae. aegypti"])
