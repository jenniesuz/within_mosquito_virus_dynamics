library(here)
library(binom)
library(bbmle)
library(ggplot2)
library(adaptivetau)
#****model code****
source(here("StochModInfectionComparisonNormal//INFModelFunc.R"))
source(here("StochModInfectionComparisonNormal//INFRepeatModelFunc.R"))
#***likelihood****
source(here("StochModInfectionComparisonNormal//INFfittingFuncsNLLTEST.R"))

#************************parameters****************************
virus_params <- function(   muV = 0.1
                            ,infRate1 = 10^-8
                            ,prodRate1 = 10            # note these parameters don't matter for infection
                            ,cellSpread1 = 10^-4        
                            ,escapeRate1 = 0.05         #
                            ,cMax = 400                 #
)
  return(as.list(environment()))
#*****************************************************************************

library(plyr)
#****generate data using 'true' parameter values****
# 30 mosquitoes, 3 virus concentrations
testDatAeg <- repeatInfModel(x=virus_params(infRate1=10^-8.5)
                          ,startingVirus=c(10^5,10^6,10^7)
)
names(testDatAeg) <- c("NumInf","ITotal","ConcMax")


testDatAlb <- repeatInfModel(x=virus_params(infRate1=10^-7.5)
                             ,startingVirus=c(10^5,10^6,10^7)
)
names(testDatAlb) <- c("NumInf","ITotal","ConcMax")

testDatAeg$Moz <- "Ae. aegypti"
testDatAlb$Moz <- "Ae. albopictus"

aegCI <- binom.confint(testDatAeg$NumInf,testDatAeg$ITotal,method="exact")
albCI <- binom.confint(testDatAlb$NumInf,testDatAlb$ITotal,method="exact")
testDatAeg$mean <- aegCI$mean
testDatAeg$lower <- aegCI$lower
testDatAeg$upper <- aegCI$upper
testDatAlb$mean <- albCI$mean
testDatAlb$lower <- albCI$lower
testDatAlb$upper <- albCI$upper
testDat <- rbind.data.frame(testDatAeg,testDatAlb)

#****plot simulated data***
ggplot(testDat) +
  geom_point(aes(x=ConcMax,y=mean,col=Moz)) +
  geom_errorbar(aes(x=ConcMax,ymin=lower,ymax=upper,col=Moz))


#*****************************************************************************
virus_params <- function(   muV = 0.1
                            ,infRate1 = 10^-8.5
                            ,infRate2 = 10^-8
                            ,prodRate1 = 1              # note these parameters don't matter for infection
                            ,cellSpread1 = 10^-4        
                            ,escapeRate1 = 0.05         #
                            ,cMax = 400                 #
)
  return(as.list(environment()))
#******************************************************************

start <- Sys.time()
trace<-3
#********initial parameter values****
init.pars.fit <- c(
  log_infRate1=log(10^-8)
  ,log_infRate2=log(10^-8)
)

#********optimise*******
optim.vals <- optim(par = init.pars.fit
                    , objFXN
                    , fixed.params = virus_params()
                    , dat = testDat
                    , nSimulations = 30
                    , control = list(trace = trace
                                     #,abstol=0.05
                                     #,reltol=0.05
                                     ,maxit=200
                    )
                    , method ="SANN" #"Nelder-Mead" #"SANN" 
)

end <- Sys.time()
end-start

# Done using tau leap epilson 0.005. Takes c. 24hrs for first SANN

# init params: 10^-8, 10^-8
# SANN 1 diff parms: 
# $par
#log_infRate1 log_infRate2 
#-19.66394    -17.57241 
# $value [1] initial 41.789714 final 6.293457

# SANN 2 diff parms:
#initial       value 7.111196
#final         value 5.393974
#$par
#log_infRate1 log_infRate2 
#-19.97774    -17.24021 

# Nelder Mead, abstol 0.05, reltol = 0.05, didn't converge estimated parameters: 
#log_infRate1 log_infRate2 
#-20.08099    -17.36464 

# repeated and converged
#log_infRate1 log_infRate2 
#-19.86871    -17.44994 
#$value
#[1] 6.337965

# very similar to initial SANN - just run with SANN?

# try with epilon 0.05 and see time reduction - only takes 18 mins
# sann objective function values
#initial       value 60.396833
#final         value 5.784653
#log_infRate1 log_infRate2 
#-8.673151    -7.548753 
