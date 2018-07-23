#' Run the SIR algorithm for a draws and stock (temporary)
#'
#' @param draws data.frame of prior draws for SIR
#' @param deplete.mean Final year depletion (not in log space)
#' @param deplete.cv Final year CV, used as SD
#' @param deplete.distribution The likelhiood to use, either (1) normal or
#'   (2) lognormal
#' @param pct.keep The percentage of "keepers" from the total. Default 10\%.
#' @param harvest.sd The mean and SD of the terminal year
#'   penalty on fishing pressure. If either is NULL it is ignored.
#' @param harvest.mean Same as harvest.sd but the mean.
#' @param harvest.distribution Same choice as deplete.distribution.
#' @param Catch Vector of catches, one for each year.
#' @param ProcessError Flag for whether to include process error
#'   deviations, passed on to model. Defaults to TRUE.
#' @param penalties A list specifying the penalties to use.
#' @param simulation See documentation for AgeModel.
#' @return A list containing depletion, SSB, and harvest rate (U) for
#'   posterior draws, and a vector of Keepers
#'
#' @export
#'
run.SIR <- function(Catch, draws, deplete.mean=NULL, deplete.cv=NULL,
                    deplete.distribution=1, harvest.distribution=1,
                    harvest.mean=NULL, harvest.sd=NULL, pct.keep=10,
                    ProcessError=TRUE, penalties=NULL, simulation=NULL){
  ## store results
  nrep <- nrow(draws)
  NY <- length(Catch)
  Bstore <- array(dim=c(NY,nrep))
  Dstore <- array(dim=c(NY,nrep))
  Ustore <- array(dim=c(NY,nrep))
  Bscaledstore <- Uscaledstore <- array(dim=c(NY,nrep))
  LikeStore <- rep(0, len=nrep)
  bmsy <- umsy <- cmsy <- crashed <- rep(NA, len=nrep)
  B <- array(dim=(NY+1))  #stock biomass
  ## Recruitment deviations (rows are replicates, columns years)
  recdevs <- matrix(NA, nrow=nrep, ncol=NY)
  ## #########################
  ## run iterations
  ## #########################
  draws <- as.matrix(draws) # less overhead than data.frame for subsetting
  for (irep in 1:nrep){   ##loop over replicates
    ## set priors
    InitialDeplete <- draws[irep, 'InitialPrior']
    Carry <- draws[irep,'Cprior']
    ## life history
    Steep <- draws[irep,'h']
    NatMort <- draws[irep,'M']
    ## Sometimes these can round down to 0 so arbitrarily setting a min
    AgeMat <- max(1,round(draws[irep,'tm']))
    AgeMax <- max(4, AgeMat,ceiling(draws[irep,'tmax']))
    Sigma <- draws[irep,'Sigma']
    ## growth
    lwa <- draws[irep,'lwa']
    lwb <- draws[irep,'lwb']
    Linf <- draws[irep,'Loo']
    vbk <- draws[irep,'K']
    ages <- 1:AgeMax
    Length <- Linf*(1-exp(-vbk*ages))
    Weight <- lwa * Length ^ lwb
    out <- AgeModel(Catch, AgeMat, Steep,NatMort, AgeMax, Carry,
                    Weight,InitialDeplete,Sigma, ProcessError=ProcessError,
                    simulation=simulation)
    pop <- out$pop
    ## plot(Years,pop,type="l",ylim=c(0,max(pop)))
    Deplete <- pop/Carry
    ## Calculate U/UMSY and B/BMSY by year
    uscaled <- out$hr/out$umsy
    bscaled <- out$Vpop/out$bmsy
    ## This is the final year "penalty" on an assumed level of depletion or
    ## harvest rate (or both). If the biomass crashes it will return a
    ## vector with NA's, so catch those here and assign a zero likelihood
    if(all(!is.na(pop))){
      loglike <- 0
      if(!is.null(deplete.mean) & !is.null(deplete.cv)){
        x <- ifelse(deplete.distribution==1, bscaled[NY], log(bscaled[NY]))
        loglike <- loglike + dnorm(x,mean=deplete.mean,sd=deplete.cv, log=TRUE)
      }
      ## If specified, include penalty for harvest rate in final year
      if(!is.null(harvest.mean) & !is.null(harvest.sd)){
        x <- ifelse(harvest.distribution==1, uscaled[NY], log(uscaled[NY]))
        loglike <- loglike + dnorm(x, mean=harvest.mean, sd=harvest.sd, log=TRUE)
      }
      LikeStore[irep] <- exp(loglike)
    } # otherwise it is 0 by default
    Bstore[,irep] <- pop
    Dstore[,irep] <- Deplete
    Ustore[,irep] <- out$hr
    Uscaledstore[,irep] <- uscaled
    Bscaledstore[,irep] <- bscaled
    cmsy[irep] <- out$cmsy
    umsy[irep] <- out$umsy
    bmsy[irep] <- out$bmsy
    crashed[irep] <- out$crashed
    recdevs[irep,] <- out$devs
  } #end of loop over replicates

  ## filter keepers
  Nkeep <- floor((pct.keep/100*nrep))
  CumLike <- c(0,cumsum(LikeStore))
  BreakPoints <- runif(Nkeep,min=0,max=CumLike[nrep+1])
  get.keepers <- function(Nkeep, CumLike, BreakPoints){
    ## This function replaces Ray's old code with a more efficienct version
    ## that loops over the BreakPoints instead of the reps. Probably a
    ## better way to do this but OK for now.  It returns a vector of
    ## indices for which draws were "kept"
    Keepers <- rep(NA, len=Nkeep)
    for(i in 1:Nkeep){
      Keepers[i] <- which(CumLike > BreakPoints[i])[1]-1
    }
    return(sort(Keepers))
  }
  ## Keepers <- array(dim=(Nkeep))
  ## k <- 1
  ## for (i in 1:nrep){  #find the keepers
  ##   Lower <- CumLike[i]; Upper <- CumLike[i+1]
  ##   j <- which(BreakPoints > Lower & BreakPoints < Upper)
  ##   Nhits <- length(j)
  ##   if (Nhits>0) {Keepers[k:(k+Nhits-1)] <-  i; k <- k+Nhits}
  ## }
  K <- get.keepers(Nkeep, CumLike, BreakPoints)
  print(paste0("# crashed=", sum(crashed),";  # unique draws=",length(unique(K))))
  penalties <- list(deplete.mean=deplete.mean,
                deplete.cv=deplete.cv, deplete.distribution,
                harvest.mean=harvest.mean,
                harvest.sd=harvest.sd, harvest.distribution=harvest.distribution)
  return(list(ssb=Bstore[,K], U=Ustore[,K], depletion=Dstore[,K],
              umsy=umsy[K], cmsy=cmsy[K], bmsy=bmsy[K], uscaled=Uscaledstore[,K],
              bscaled=Bscaledstore[,K], recdevs=recdevs, penalties=penalties,
              Keepers=K, crashed=crashed, likes=LikeStore, Catch=Catch,
              draws=draws))
}
