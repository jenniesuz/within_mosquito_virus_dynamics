library(mvtnorm)
library(gridExtra)
library(ggplot2)
library(binom)
library(beepr)
library(parallel)
library(here)
library(bbmle)
#****model code****
source(here("StochModDisseminationComparisonNormal//InfDissModelFuncs.R"))
source(here("StochModDisseminationComparisonNormal//InfDissRepeatModelFunc.R"))
#***likelihood****
source(here("StochModDisseminationComparisonNormal//DISSfittingFuncsNLL.R"))

#******************data to fit to********************
competenceDat <- read.csv(here(".//Data//datCiotaOnyango.csv"))
competenceDat <- competenceDat[competenceDat$Ref %in% "Onyango 2020",]
bin <- binom.confint(x=competenceDat$NumDiss,n=competenceDat$NumInf,method="exact")
competenceDat$meanDiss <- bin$mean
competenceDat$lowerDiss <- bin$lower
competenceDat$upperDiss <- bin$upper
#****************************************************

#*******************parameters*****************
virus_params <- function(   muV = 0.1
                            ,infRate = 10^-9 
                            ,prodRate = 170  
                            ,cellSpread = 0.00023 
                            ,escapeRate = 0.25 
                            ,cMax = 400  
                            ,hMax= 900
)
return(as.list(environment()))



dat <- competenceDat[competenceDat$Moz %in% "Ae. albopictus",]

fit <- readRDS("DissModelFitSepAllParms230608.rds")

nll.binom.mle2()

#saveRDS(llsdp,"albLLRange.230603")

# log_infRate1 log_infRate2 
# -20.50493    -16.78952   


fitMLE <- mle2(nll.binom.mle2,start=list(prodRate=fit$par[2],cellSpread=fit$par[4],escapeRate=fit$par[6])
               ,fixed=list(muV=0.1,infRate=exp(-16.78952)))

