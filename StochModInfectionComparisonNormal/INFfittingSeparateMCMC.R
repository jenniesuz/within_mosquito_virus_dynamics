library(here)
library(binom)
library(ggplot2)
library(bbmle)
library(parallel)
library(adaptivetau)
library(mvtnorm)
#****model code****
source(here("StochModInfectionComparisonNormal//INFmodelFunc.R"))

source(here("StochModInfectionComparisonNormal//INFRepeatModelFunc.R"))
#***likelihood****
source(here("StochModInfectionComparisonNormal//INFfittingFuncsNLL.R"))

require(boot); require(deSolve); require(ellipse); require(coda); require(parallel); require(mnormt); require(emdbook)

#******************data to fit to********************
competenceDat <- read.csv(here("StochModInfectionComparisonNormal//datCiotaOnyango.csv"))
competenceDat <- competenceDat[competenceDat$Ref %in% "Ciota 2017",]
bin <- binom.confint(x=competenceDat$NumInf,n=competenceDat$ITotal,method="exact")
competenceDat$meanInf <- bin$mean
competenceDat$lowerInf <- bin$lower
competenceDat$upperInf <- bin$upper
#****************************************************


## Function that makes a list of disease parameters with default values
virus_params <- function(   muV = 0.1
                            ,infRate1 = 10^-8
                            ,infRate2 = 10^-8
                            ,prodRate1 = 10            # note these parameters don't matter for infection
                            ,cellSpread1 = 10^-4        
                            ,escapeRate1 = 0.05         #
                            ,cMax = 400                 #
)
  return(as.list(environment()))


## Log-Prior (assume uninformative)
lprior <- function(parms=virus_params()) with(parms, {
    lp <- 0     ## whatever the parameters are, assume they have the same probability
    return(lp)
})

## Convenience function that sums log-likelihood & log-prior for
## evaluation inside MCMC sampler.
llikePrior <- function(fit.params=NULL, ## parameters to fit
                       ref.params = virus_params(), ## reference parameters
                       dat=competenceDat) { ## observed data
    parms <- within(ref.params, { ## subs fitting parameters into reference parameter vector
        for(nm in names(fit.params)) assign(nm, as.numeric(fit.params[nm]))
        rm(nm)
    })
    -nll.binom(parms, dat=dat,nSims=30) + lprior(parms)
}
llikePrior(dat=competenceDat)




## Want to be able to easily log and unlog parameters
logParms <- function(fit.params) {
    fit.params <- log(fit.params)
    names(fit.params) <- paste0('log',names(fit.params))
    return(fit.params)
}
logParms(c(infRate1 = 10^-8, infRate2=10^-8))

unlogParms <- function(fit.params) {
    fit.params <- exp(fit.params)
    names(fit.params) <- sub('log','', names(fit.params))
    return(fit.params)
}
unlogParms(logParms(c(infRate1=10^-8, infRate2=10^-8)))

## set bounds on initial parameter guesses
initBounds <- data.frame(rbind( ## for initial conditions
                               c(log(10^-10),log(10^-6)) ## infRate1
                               ,c(log(10^-10), log(10^-6)) ## infRate2
                               ,c(log(0.01),log(0.3)))) ## muV
colnames(initBounds) <- c('lower','upper')
rownames(initBounds) <- c('loginfRate1','loginfRate2','logmuV')
class(initBounds[,2]) <- class(initBounds[,1]) <- 'numeric'
initBounds

##  randomly select a value that is uniformly distributed between these bounds
initRand <- function(fit.params) {
    fit.params <- logParms(fit.params)
    tempnm <- names(fit.params)
    for(nm in tempnm) fit.params[nm] <- runif(1, min = initBounds[rownames(initBounds)==nm, 'lower'], 
                                              max =  initBounds[row.names(initBounds)==nm, 'upper'])
    return(unlogParms(fit.params))
}
## give it parameter vector and it will simulate random values within the bounds for those values
initRand(c(infRate1=10^-8, infRate2=10^-8,muV=0.1)) 

## Flexible Metropolis-Hastings Sampler
mcmcSampler <- function(init.params, ## initial parameter guess
                        randInit = T, ## if T then randomly sample initial parameters instead of above value
                        seed = 1, ## RNG seed
                        ref.params=virus_params(), ## fixed parameters
                        obsDat = competenceDat, ## data
                        proposer = sequential.proposer(sdProps=sdProps), ## proposal distribution
                        niter = 100, ## MCMC iterations
                        nburn = 0, ## iterations to automatically burn
                        adaptiveMCMC = F, ## adapt proposal distribution?
                        startAdapt = 150, ## start adapting at what iteration?
                        adptBurn = 200, ## ignore first so many iterations for adapting posterior
                        verbose=0, ## if >2 browses, if >1 prints progress
                        tell = 100) { ## how often to print progress
    if(verbose>2) browser()
    if(randInit) init.params <- initRand(init.params)
    current.params <- init.params
    nfitted <- length(current.params) ## number fitted parameters
    vv <- 2 ## mcmc iteration (started at 1 so we're already on 2
    accept <- 0 ## initialize proportion of iterations accepted
    ## Calculate log(likelihood X prior) for first value
    curVal <- llikePrior(current.params, ref.params = ref.params, dat=obsDat)
    ## Initialize matrix to store MCMC chain
    out <- matrix(NA, nr = niter, nc=length(current.params)+1)
    out[1,] <- c(current.params, ll = -curVal) ## add first value
    colnames(out) <- c(names(current.params), 'll') ## name columns
    ## Store original covariance matrix
    if(proposer$type=='block') originalCovar <- get('covar', envir = environment(proposer$fxn)) 
    while(vv <= niter) {
        if ((verbose > 1) || (verbose && (vv%%tell == 0))) print(paste("on iteration",vv,"of", niter + 1))
        ## Adaptive MCMC: adapt covariance every 50 iterations (don't
        ## do it more often because it adds to coputational burden.
        if(adaptiveMCMC & proposer$type=='block' & vv > startAdapt & vv %% 50 == 0) {
            adptBurn <- min((startAdapt-50), adptBurn)
            ## Below equation gives ideal covariance-variance matrix based on posterior
            adaptedCovar <- 2.38^2 / nfitted * cov.wt(log(out[adptBurn:(vv-1),1:nfitted]))$cov
            ## Take a weighted average of the original & the empirical cov-var matrices to ensure
            ## that we never let the matrix collapse to zero (ie if the empirical one is zero
            ## because we haven't accepted anything yet)
            adaptedCovar <- adaptedCovar*.95 + originalCovar*.05 ## 95% adapted & 5% original
            rownames(adaptedCovar) <- colnames(adaptedCovar) <- names(current.params)
            assign('covar', adaptedCovar, envir = environment(proposer$fxn))
        }
        proposal <- proposer$fxn(logParms(current.params))
        proposal <- unlogParms(proposal)
        propVal <- llikePrior(proposal, ref.params = ref.params, dat=obsDat)
        lmh <- propVal - curVal ## likelihood ratio = log likelihood difference
        if (is.na(lmh)) { ## if NA, print informative info but don't accept it
            print(list(lmh=lmh, proposal=exp(proposal), vv=vv, seed=seed))
        } else { ## if it's not NA then do acception/rejection algorithm
            if (verbose > 1) print( c(lmh=lmh, propVal=propVal) )
            ## if MHR >= 1 or a uniform random # in [0,1] is <= MHR, accept otherwise reject
            if ( (lmh >= 0) | (runif(1,0,1) <= exp(lmh)) ) {
                current.params <- proposal
                if (vv>nburn) accept <- accept + 1 ## only track acceptance after burn-in
                curVal <- propVal
            }
        }
        out[vv, ] <- c(current.params, ll=curVal)
        vv <- vv+1
        aratio <- accept/((vv-nburn))
    }
    colnames(out) <- c(names(current.params), 'll')
    samp <- as.mcmc(out[1:nrow(out)>(nburn+1),])
    return(list(ref.params=ref.params
              , seed = seed
              , init.params = init.params
              , aratio = aratio
              , samp = samp
                ))
}

## Sequential proposal function: Propose one parameter at a time
sequential.proposer <- function(sdProps) {
    nfitted <- length(sdProps)
    on <- 0
    return(list(sdProps = sdProps, type = 'sequential',
                fxn = function(current) {
                    proposal <- current
                    proposal[on + 1] <- proposal[on + 1] + rnorm(1, mean = 0, sd = sdProps[on + 1])
                    on <<- (on+1) %% nfitted
                    proposal
                }))
}


## Propose parameters within blocks
multiv.proposer <- function(covar, blockLS = list(rownames(covar))) {
    nblocks <- length(blockLS)
    on <- 0
    return(list(type = 'block',
                fxn = function(current) {
                    proposal <- current + rmnorm(1, mean = 0, varcov = covar)
                    propsosal <- as.vector(proposal)
                    names(proposal) <- names(current)
                    proposal
                }))
}







samp_Seq <- mcmcSampler(init.params = c(infRate1=10^-9,infRate2=10^-7,muV=0.1)
                      , seed = 1
                      , proposer = sequential.proposer(sdProps=c(.15,.15))
                      , randInit = T
                      , niter = 500)


class(samp_Seq$samp)
## The coda package already knows how to plot MCMC objects by default
par('ps'=4, las = 0.2)
plot(samp_Seq$samp[,2])













####################################################################################################
## Adaptive proposals
####################################################################################################
mcmcParams <- within(mcmcParams, {
    niter <- 2000 ## let's increase the # of iterations
    verbose <- 1
    tell <- 100
})
mcmcParams_Adaptive <- within(mcmcParams, {
                      proposer <- multiv.proposer(covar=matrix(c(.1,0,0,.1),2,2))
                  })

## run4 <- doChains(1:4, mcmcParams) ## do two chains with seeds at 1:2
## run4A <- doChains(1:4, mcmcParams_Adaptive) ## do two chains with seeds at 1:2
## save(run4, run4A, file = 'MCMC_SI_runs.Rdata')
load(file = 'MCMC_SI_runs.Rdata')

run4$aratio
run4A$aratio

par(oma = c(0,0,2,0), bty='n', 'ps' = 18)
plot(run4$chains)
mtext('Sequential Sampling', side = 3, outer = T, line = 0)
 
plot(run4A$chains)
mtext('Adaptive Sampling', side = 3, outer = T, line = 0)

graphics.off()

gelman.diag(run4$chains[,c('alpha','Beta')])
gelman.diag(run4A$chains[,c('alpha','Beta')])

summary(run4$chains) ## Posterior credible intervals
summary(run4A$chains)

par(mar = c(5,6,1,1), las = 1, 'ps' = 18, mfrow = c(2,1))
for(nm in c('4','4A')) {
    res <- get(paste0('run', nm))$chains
    plot(unlist(res[,'alpha']), unlist(res[,'Beta']),
         xlab = expression(alpha),
         ylab = expression(beta),
         log = 'xy',
         type = 'p',
         cex = .7, pch = 16,
         col = gray(.5, alpha = .1),
         xlim = c(2,15),
         ylim = c(.3, 2))
    ## Bayesian 95% credible contour calculated by finding highest posterior density region.
    HPDregionplot(res, 1:2,
                  prob = .95,
                  n = 40,
                  lwd = 2,
                  col = 'red',
                  add = T) ## add to current plot
}

