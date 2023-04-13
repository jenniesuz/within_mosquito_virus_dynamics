library(adaptivetau)

infectionModel <- function(startingVirus
                           ,muV = 0.1
                           ,infRate = 10^-7.5
                           ,prodRate =   1    
                           ,cellSpread = 10^-4  
                           ,escapeRate = 0.05
                           ,cMax = 400 
){
  
  
  params <- list(muV = muV
                 ,probInf = infRate
                 ,prodRate =   prodRate 
                 ,cellSpread = cellSpread
                 ,escapeRate = escapeRate    
                 ,cMax = cMax            
  )
  
  transitions <- list(c(Gv = -1,Mc = +1)
                      ,c(Gv = -1)
                      ,c(Mc = +1)
                      ,c(Mv = +1)
                      ,c(Mv = -1)
                      ,c(Mv = -1)
  )
  
  
  lvrates <- function(y,params,t){
    return( c(y["Gv"]*params$probInf*(params$cMax-y["Mc"])
              ,params$muV*y["Gv"]
              ,params$cellSpread*y["Mc"]*(params$cMax-y["Mc"])
              ,params$prodRate*y["Mc"]
              ,params$muV*y["Mv"]
              ,params$escapeRate*y["Mv"]
    )
    
    )
  }
  
  out<-ssa.adaptivetau(c(Gv = round(startingVirus*0.003,0), Mc = 0, Mv = 0),
                       transitions, lvrates, params, tf=120
                       , tl.params=list(epsilon=0.005)) 
  return(data.frame(out))
}



