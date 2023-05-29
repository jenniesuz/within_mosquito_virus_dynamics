library(bbmle)
dat=competenceDat[competenceDat$Moz %in% "Ae. aegypti",]

tsum <- summary(test)

repFits <- lapply(1:2,function(x){
 rep <- x
  test <- mle2(function(logmuV,loginfRate){nll.binom.mle2()}             # fit first model
     ,start=list(logmuV=log(0.1),loginfRate=log(10^-8.905177))
     ,skip.hessian=FALSE,hessian.opts=c(method="complex"))  
  return(c(x,coef(test)))
})
do.call(rbind,repFits)

