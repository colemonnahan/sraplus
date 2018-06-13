
#' Run the SIR algorithm for a draws and stock (temporary)
#' @param draws data.frame of prior draws for SIR
#' @param deplete.mean Final year depletion (not in log space)
#' @param deplete.cv Final year CV, used as SD
#' @param pct.keep The percentage of "keepers" from the total. Default 10%.
#' @param harvest.sd, harvest.mean The mean and SD of the terminal year
#'   penalty on fishing pressure. If either is NULL it is ignored.
#' @param Catch Vector of catches, one for each year.
#' @return A list containing depletion, SSB, and harvest rate (U) for
#'   posterior draws, and a vector of Keepers
run.SIR <- function(Catch, draws, deplete.mean=NULL, deplete.cv=NULL,
                    harvest.mean=NULL, harvest.sd=NULL, pct.keep=10){
  ## store results
  nrep <- nrow(draws)
  NY <- length(Catch)
  Bstore <- array(dim=c(NY,nrep))
  Dstore <- array(dim=c(NY,nrep))
  Ustore <- array(dim=c(NY,nrep))
  LikeStore <- array(dim=nrep)
  B <- array(dim=(NY+1))  #stock biomass
  ## #########################
  ## run iterations
  ## #########################
  for (irep in 1:nrep){   ##loop over replicates
    ## set priors
    InitialDeplete <- draws$InitialPrior[irep]
    Carry <- draws$Cprior[irep]
    ## life history
    Steep <- draws$h[irep]
    NatMort <- draws$M[irep]
    ## Sometimes these can round down to 0 so arbitrarily setting a min
    AgeMat <- max(1,round(draws$tm[irep],0))
    AgeMax <- max(2, AgeMat,ceiling(draws$tmax[irep]))
    Sigma <- draws$Sigma[irep]
    ## growth
    lwa <- draws$lwa[irep]
    lwb <- draws$lwb[irep]
    Linf <- draws$Loo[irep]
    vbk <- draws$K[irep]
    ages <- 1:AgeMax
    Length <- Linf*(1-exp(-vbk*ages))
    Weight <- lwa * Length ^ lwb
    out <- AgeModel(Catch, AgeMat, Steep,NatMort, AgeMax, Carry, Weight,InitialDeplete,Sigma)
    pop <- out$pop
    ## plot(Years,pop,type="l",ylim=c(0,max(pop)))
    Deplete <- pop/Carry
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
        loglike2 <- dnorm(out$hr[NY],mean=harvest.mean,sd=harvest.sd, log=TRUE)
      like <- exp(loglike1+loglike2)
    }
    Bstore[,irep] <- pop
    Dstore[,irep] <- Deplete
    Ustore[,irep] <- out$hr
    LikeStore[irep] <- like
  } #end of loop over replicates

  ## filter keepers
  Nkeep <- floor((pct.keep/100*nrep))
  CumLike <- c(0,cumsum(LikeStore))
  BreakPoints <- runif(Nkeep,min=0,max=CumLike[nrep+1])
  Keepers <- array(dim=(Nkeep))
  k <- 1
  for (i in 1:nrep){  #find the keepers
    Lower <- CumLike[i]; Upper <- CumLike[i+1]
    j <- which(BreakPoints > Lower & BreakPoints < Upper)
    Nhits <- length(j)
    if (Nhits>0) {Keepers[k:(k+Nhits-1)] <-  i; k <- k+Nhits}
  }
  print(paste("% unique draws=",100*length(unique(Keepers))/Nkeep))

  return(list(depletion=Dstore[, Keepers], ssb=Bstore[,Keepers],
              U=Ustore[,Keepers], final=list(deplete.mean=deplete.mean,
              deplete.cv=deplete.cv, harvest.mean=harvest.mean,
              harvest.sd=harvest.sd), Keepers=Keepers))
}
