#' Run the SIR algorithm for a draws and stock (temporary)
#'
#' @param draws data.frame of prior draws for SIR
#' @param deplete.mean Final year depletion (not in log space)
#' @param deplete.cv Final year CV, used as SD
#' @param pct.keep The percentage of "keepers" from the total. Default 10\%.
#' @param harvest.sd The mean and SD of the terminal year
#'   penalty on fishing pressure. If either is NULL it is ignored.
#' @param harvest.mean Same as harvest.sd but the mean.
#' @param Catch Vector of catches, one for each year.
#' @param ProcessError Flag for whether to include process error
#'   deviations, passed on to model. Defaults to TRUE.
#' @return A list containing depletion, SSB, and harvest rate (U) for
#'   posterior draws, and a vector of Keepers
#'   @export
run.SIR <- function(Catch, draws, deplete.mean=NULL, deplete.cv=NULL,
                    harvest.mean=NULL, harvest.sd=NULL, pct.keep=10,
                    ProcessError=TRUE, simulation=NULL){
  ## store results
  nrep <- nrow(draws)
  NY <- length(Catch)
  Bstore <- array(dim=c(NY,nrep))
  Dstore <- array(dim=c(NY,nrep))
  Ustore <- array(dim=c(NY,nrep))
  Uscaledstore <- array(dim=c(NY,nrep))
  LikeStore <- array(dim=nrep)
  umsy <- cmsy <- rep(NA, len=nrep)
  B <- array(dim=(NY+1))  #stock biomass
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
    AgeMax <- max(2, AgeMat,ceiling(draws[irep,'tmax']))
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
    ## Calculate U/UMSY by year
    uscaled <- out$hr/out$umsy
    ## This is the final year "penalty" on an assumed level of depletion or
    ## harvest rate (or both). If the biomass crashes it will return a
    ## vector with NA's, so catch those here and assign a zero likelihood
    if(any(is.na(pop))){
      like <- 0
    } else {
      loglike1 <- loglike2 <- 0
      if(!is.null(deplete.mean) & !is.null(deplete.cv))
        loglike1 <- dnorm(Deplete[NY],mean=deplete.mean,sd=deplete.cv, log=TRUE)
      ## If specified, include penalty for harvest rate in final year
      if(!is.null(harvest.mean) & !is.null(harvest.sd))
        loglike2 <- dnorm(uscaled[NY],mean=harvest.mean,sd=harvest.sd, log=TRUE)
      like <- exp(loglike1+loglike2)
    }
    Bstore[,irep] <- pop
    Dstore[,irep] <- Deplete
    Ustore[,irep] <- out$hr
    Uscaledstore[,irep] <- uscaled
    LikeStore[irep] <- like
    cmsy[irep] <- out$cmsy
    umsy[irep] <- out$umsy
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
  print(paste("% unique draws=",100*length(unique(K))/Nkeep))
  final <- list(deplete.mean=deplete.mean,
                deplete.cv=deplete.cv, harvest.mean=harvest.mean,
                harvest.sd=harvest.sd)
  return(list(depletion=Dstore[, K], ssb=Bstore[,K],
              U=Ustore[,K], Keepers=K, final=final, umsy=umsy[K],
              cmsy=cmsy[K], uscaled=Uscaledstore[,K]))
}
