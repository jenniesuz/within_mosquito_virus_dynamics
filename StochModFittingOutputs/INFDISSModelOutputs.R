library(gridExtra)
library(ggplot2)
library(binom)
library(viridis)

aegInf <- readRDS("aegInffits.rds") 
albInf <- readRDS("albInffits.rds")

aegInfParams <- sapply(aegInf,"[[",2)
albInfParams <- sapply(albInf,"[[",2)

aegInfParams <- cbind.data.frame(muV=aegInfParams[1,]        # convert from list to dataframe
                           ,probInf=aegInfParams[2,]
                           ,species="Ae. aegypti"
)

aegInfParams$run <- 1:100


aegInfFit <- sapply(aegInf,"[[",1)
aegInfHess <- aegInfFit[6,]
aegInfLik <- aegInfFit[2,]
aegInfLik <- do.call(rbind.data.frame,aegInfLik)
names(aegInfLik) <- "Likelihood"
aegInfLik$run <- 1:100


albInfFit <- sapply(albInf,"[[",1)


#**********Confidence intervals*********
#*#*****************confidence intervals for parameter estimates*********
confFromHessFunc <- function(hessianMatrix,parms){

  fisherInfMatrix <- solve(hessianMatrix) ## invert the Hessian, to estimate the covar-var matrix of parameter estimates

  # Finds the critical z value
  conf.level <- 0.95
  crit <- qnorm((1 + conf.level)/2)

  ci1 <- log(as.numeric(parms[1])) + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[1, 1]))

  ci2 <- log(as.numeric(parms[2])) + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[2, 2]))

return(rbind(exp(ci1),exp(ci2)))
}


confFromHessFunc(aegInfHess[[1]],parms=(aegInfParams[1,]))



#***********************************

# read in model fits fits
albF <- readRDS("albDissFits220418.rds")
aegF <- readRDS("aegDissFits220418.rds")

# function to change to dataframe and then bind together
listToDatParms <- function(output=alb, species="alb"){
  params <- lapply(1:length(output),function(x){
    temp <- output[[x]]
    if(is.list(temp)==T){
      temp <- temp[[1]][1]
      temp <- do.call(c,temp)
      names(temp) <- c("prodRate","cellSpread","escapeRate")
      return(exp(temp))
    }else{
      temp <- c(NA,NA,NA)
      names(temp) <- c("prodRate","cellSpread","escapeRate")
      return(temp)
    }
  })
  params <- do.call(rbind.data.frame,params)
  names(params) <- c("prodRate","cellSpread","escapeRate")
  params$species=species
  return(params)
}

albDissParams <- listToDatParms(albF)
aegDissParams <- listToDatParms(output=aegF,species="aeg")

# function to change to dataframe and then bind together
listToDatHess<- function(output=alb, species="alb"){
  hessians <- lapply(1:length(output),function(x){
    temp <- output[[x]]
    if(is.list(temp)==T){
      hess <- temp[[1]][6]
      return(hess)
    }else{
      hess <- c(NA)
      return(temp)
    }
  })
  return(hessians)
}

albDissHess <- listToDatHess(albF)
aegDissHess <- listToDatHess(aegF)


#**********Confidence intervals*********
#*#*****************confidence intervals for parameter estimates*********
confFromHessFunc <- function(hessianMatrix,parms){
  
  fisherInfMatrix <- solve(hessianMatrix) ## invert the Hessian, to estimate the covar-var matrix of parameter estimates
  
  # Finds the critical z value
  conf.level <- 0.95
  crit <- qnorm((1 + conf.level)/2)
  
  ci1 <- log(as.numeric(parms[1])) + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[1, 1]))
  
  ci2 <- log(as.numeric(parms[2])) + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[2, 2]))
  
  ci3 <- log(as.numeric(parms[3])) + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[3, 3]))
  
  
  return(rbind(exp(ci1),exp(ci2),exp(ci3)))
}


confFromHessFunc(albDissHess[[1]]$hessian,parms=(albDissParams[1,]))



#***********************************
