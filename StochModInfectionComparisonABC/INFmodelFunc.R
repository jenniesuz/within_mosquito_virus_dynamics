library(adaptivetau)

infectionModel <- function(startingVirus,parms=c(par1=0.03,par2=10^-8)
                     ){
  
  
  params <- list(muV = as.numeric(parms[1])
                 ,probInf = as.numeric(parms[2])
                 ,prodRate =   1    
                 ,cellSpread = 10^-4  
                 ,escapeRate = 0.05    
                 ,cMax = 400            
  )
  
  transitions <- list(c(Gv = -1)
                      ,c(Gv = -1,Mc = +1)
                      ,c(Mc = +1)
                      ,c(Mv = +1)
                      ,c(Mv = -1)
                      ,c(Mv = -1)
  )
  
  
  lvrates <- function(x,params,t){
    return( c(params$muV*x["Gv"]
              ,x["Gv"]*params$probInf*(params$cMax-x["Mc"])
              ,params$cellSpread*x["Mc"]*(params$cMax-x["Mc"])
              ,params$prodRate*x["Mc"]
              ,params$muV*x["Mv"]
              ,params$escapeRate*x["Mv"]
    )
    
    )
  }
  
    out<-ssa.adaptivetau(c(Gv = round(startingVirus*0.003,0), Mc = 0, Mv = 0),
                      transitions, lvrates, params, tf=120) 
    return(data.frame(out))
 }

