library(gridExtra)
library(ggplot2)
library(binom)
library(beepr)

source("3_INFfittingFuncsNLL.R")
source("3_INFmodelFunc.R")

virus_params <- function(   muV = 0.1
                            ,infRate = 10^-7.5
                            ,prodRate =   1    
                            ,cellSpread = 10^-4  
                            ,escapeRate = 0.05    
                            ,cMax = 400  
)
return(as.list(environment()))

test <- sim.func(10^5,parms=virus_params())

competenceDat <- read.csv("datCiotaOnyango.csv")
competenceDat <- competenceDat[competenceDat$Ref %in% "Ciota 2017",]
bin <- binom.confint(x=competenceDat$NumInf,n=competenceDat$ITotal,method="exact")
competenceDat$meanInf <- bin$mean
competenceDat$lowerInf <- bin$lower
competenceDat$upperInf <- bin$upper

#****************************AEGYPTI HND****************************************
#***************************************************************************
init.pars.fit1 <- c(
  log_muV=log(0.1)
  ,log_infRate=log(10^-9.5)
)
#**************Optimise*******************************************
trace <- 3

optim.vals <- optim(par = init.pars.fit1
                    , objFXN
                    , fixed.params = virus_params()
                    , dat = competenceDat[(competenceDat$Moz %in% "Ae. aegypti"),]
                    , control = list(trace = trace, maxit = 200, reltol = 10^-7)
                    , method = "Nelder-Mead" # 
                    , hessian = T)
optim.vals # convergence 0 means algorithm converged
beep()

MLEfits <- optim.vals$par 
exp(MLEfits)

#muV  infRate 
aegHNDmuV <- exp(MLEfits)[1]
aegHNDinfRate <- exp(MLEfits)[2]
#******************************************************************************
#****************************AEGYPTI CAM****************************************
#***************************************************************************
init.pars.fit1 <- c(
  log_muV=log(0.1)
  ,log_infRate=log(10^-9.5)
)
#**************Optimise*******************************************
trace <- 3

optim.vals <- optim(par = init.pars.fit1
                    , objFXN
                    , fixed.params = virus_params()
                    , dat = competenceDat[(competenceDat$Moz %in% "Aeae")&(competenceDat$Visolate %in% "CAM"),]
                    , control = list(trace = trace, maxit = 500, reltol = 10^-7.5)
                    , method = "Nelder-Mead" # 
                    , hessian = T)
optim.vals # convergence 0 means algorithm converged
beep()

MLEfits <- optim.vals$par 
exp(MLEfits)

aegCAMmuV <- exp(MLEfits)[1]
aegCAMinfRate <- exp(MLEfits)[2]
#******************************************************************************

#****************************ALBOPICTUS HND****************************************
#***************************************************************************
init.pars.fit1 <- c(
  log_muV=log(0.1)
  ,log_infRate=log(10^-6.5)
)
#**************Optimise*******************************************
trace <- 3

optim.vals <- optim(par = init.pars.fit1
                    , objFXN
                    , fixed.params = virus_params()
                    , dat = competenceDat[(competenceDat$Moz %in% "Aea")&(competenceDat$Visolate %in% "HND"),]
                    , control = list(trace = trace, maxit = 500, reltol = 10^-7.5)
                    , method = "Nelder-Mead" # 
                    , hessian = T)
optim.vals # convergence 0 means algorithm converged
beep()

MLEfits <- optim.vals$par 
exp(MLEfits)

#muV  infRate 
albHNDmuV <- exp(MLEfits)[1]
albHNDinfRate <- exp(MLEfits)[2]
#******************************************************************************

  #****************************ALBOPICTUS CAM****************************************
  #***************************************************************************
  init.pars.fit1 <- c(
    log_muV=log(0.1)
    ,log_infRate=log(10^-6.5)
  )
#**************Optimise*******************************************
trace <- 3

optim.vals <- optim(par = init.pars.fit1
                    , objFXN
                    , fixed.params = virus_params()
                    , dat = competenceDat[(competenceDat$Moz %in% "Aea")&(competenceDat$Visolate %in% "CAM"),]
                    , control = list(trace = trace, maxit = 500, reltol = 10^-7.5)
                    , method = "Nelder-Mead" # 
                    , hessian = T)
optim.vals # convergence 0 means algorithm converged
beep()

MLEfits <- optim.vals$par 
exp(MLEfits)

#muV  infRate 
albCAMmuV <- exp(MLEfits)[1]
  albCAMinfRate <- exp(MLEfits)[2]
#******************************************************************************
  














#**************************PLOT***************************************************

#**************Aegypti***********
doseSimAegHND <- mclapply(c(10^3,10^3.5,10^4,10^4.5,10^5,10^5.5,10^6,10^6.5,10^7,10^7.5,10^8,10^8.5,10^9,10^9.5,10^10)
                          , sim.func,parms=virus_params(muV=2.436233e-01 ,infRate=1.303870e-07 ))
doseSimAegCAM <- mclapply(c(10^3,10^3.5,10^4,10^4.5,10^5,10^5.5,10^6,10^6.5,10^7,10^7.5,10^8,10^8.5,10^9,10^9.5,10^10)
                          , sim.func,parms=virus_params(muV=aegCAMmuV ,infRate=aegCAMinfRate ))


doseSimAegHND <- do.call(rbind.data.frame,doseSimAegHND)
doseSimAegCAM <- do.call(rbind.data.frame,doseSimAegCAM)

infDatAegHND <- modOutFunc(doseSimAegHND)
infDatAegCAM <- modOutFunc(doseSimAegCAM)

#**************Albo***********
doseSimAlbHND <- mclapply(c(10^3,10^3.5,10^4,10^4.5,10^5,10^5.5,10^6,10^6.5,10^7,10^7.5,10^8,10^8.5,10^9,10^9.5,10^10)
                          , sim.func,parms=virus_params(muV=albHNDmuV ,infRate=albHNDinfRate  ))
doseSimAlbCAM <- mclapply(c(10^3,10^3.5,10^4,10^4.5,10^5,10^5.5,10^6,10^6.5,10^7,10^7.5,10^8,10^8.5,10^9,10^9.5,10^10)
                          , sim.func,parms=virus_params(muV=albCAMmuV ,infRate=albCAMinfRate ))


doseSimAlbHND <- do.call(rbind.data.frame,doseSimAlbHND)
doseSimAlbCAM <- do.call(rbind.data.frame,doseSimAlbCAM)

infDatAlbHND <- modOutFunc(doseSimAlbHND)
infDatAlbCAM <- modOutFunc(doseSimAlbCAM)


aeg <- ggplot(competenceDat[(competenceDat$Moz %in% "Aeae"),]) +
  geom_point(aes(x=Conc.Min,y=meanInf,col=Visolate))+
  geom_errorbar(aes(x=Conc.Min,ymin=lowerInf,ymax=upperInf)) +
  geom_point(data=infDatAegHND,aes(x=log10(Conc.Min),y=meanInf),col="red",size=2)  +
  geom_point(data=infDatAegCAM,aes(x=log10(Conc.Min),y=meanInf),col="blue",size=2) +
  xlab("Minimum virus concentration") +
  ylab("Proportion with a midgut infection") +
 # facet_wrap(~Visolate,labeller=label_wrap_gen(width=8))+
  theme_set(theme_bw())  +    # I think this is a nicer theme for plots
  theme(panel.border = element_blank()                   # all these functions under theme can alter apperance
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.5,"line")   # here you can amend legend size and position
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,strip.text=element_text(size=5)
  )   
aeg

alb <- ggplot(competenceDat[(competenceDat$Moz %in% "Aea"),]) +
  geom_point(aes(x=Conc.Min,y=meanInf,col=Visolate))+
  geom_errorbar(aes(x=Conc.Min,ymin=lowerInf,ymax=upperInf)) +
  geom_point(data=infDatAlbHND,aes(x=log10(Conc.Min),y=meanInf),col="red",size=2)  +
  geom_point(data=infDatAlbCAM,aes(x=log10(Conc.Min),y=meanInf),col="blue",size=2) +
  xlab("Minimum virus concentration") +
  ylab("Proportion with a midgut infection") +
  # facet_wrap(~Visolate,labeller=label_wrap_gen(width=8))+
  theme_set(theme_bw())  +    # I think this is a nicer theme for plots
  theme(panel.border = element_blank()                   # all these functions under theme can alter apperance
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.5,"line")   # here you can amend legend size and position
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,strip.text=element_text(size=5)
  )   

alb
