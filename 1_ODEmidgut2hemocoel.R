library(ggplot2)
library(gridExtra)
library(deSolve)

#*********************PARAMETERS*************************************
params <- function(     
   muV =  0.1             # virus clearance/ death rate 
  ,infRate = 10^-8        # probability of infection & contact rate
  ,prodRate = 10         # virion production rate
  ,cellSpread = 10^-3.5   # rate virions in infected cells infect susceptible cells
  ,escapeRate = 0.05      # rate virions escape through basal lamina of midgut
  ,cMax = 400             # number of cells in midgut
  ,hMax = 900             # number of cells in haemocoel
 )
return(as.list(environment()))
#***************************************************************

#*****************INITIAL CONDITIONS*****************************
initial <- c(Gv = 10^6      # number of virions in bloodmeal
             ,Mci = 0       # number of infected midgut cells
             ,Mv = 0        # number of virions in midgut
             ,Hv = 0        # number of infected haemocoel cells
             ,Hci = 0        # number of virions in haemocoel
             
)

times <- seq(0,170,1)               # times to solve at
#**************************************************************

#****************MODEL*****************************************
mod <- function(tt,yy,parms) with(c(parms,as.list(yy)), {
  
  deriv <- rep(NA,5)
  
  deriv[1] <- -Gv*infRate*(cMax - Mci) - muV*Gv      # virion in bloodmeal
  
  deriv[2] <- Gv*infRate*(cMax-Mci) + cellSpread*Mci*(cMax-Mci)   # infected midgut cells
  
  deriv[3] <- prodRate*Mci - muV*Mv - escapeRate*Mv      # virions in midgut
  
  deriv[4] <- escapeRate*Mv + prodRate*Hci - muV*Hv      # virions in haemocoel
  
  deriv[5] <- Hv*infRate*(hMax-Hci)                      # infected haemocoel cells
  
  return(list(deriv))
})
#*************************************************************

#**************SIMULATE***************************************
simPop <- function(init=initial, tseq = times, modFunction=mod, parms = params()) {
  simDat <- as.data.frame(ode(init, tseq, modFunction, parms=parms))
  return(simDat)
}
#****************************************************************
