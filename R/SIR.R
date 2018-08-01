#' Run the SIR algorithm for a draws and stock (temporary)
#'
#' @param nrep Number of replicates to run.
#' @param Catch Vector of catches, one for each year.
#' @param Taxon A list with taxonomic information, including 'Class',
#'   'Order', 'Family', 'Genus', 'Species'. Any blanks will be converted to
#'   'predictive' internally when querying FishLife.
#' @param Kprior Maybe a multiplier on K?
#' @param Kscale Scalar to control the initial biomass which is generated
#' @param pct.keep The percentage of "keepers" from the total. Default
#'   10\%.
#' @param ProcessError Flag for whether to include process error
#'   deviations, passed on to model. Defaults to TRUE.
#' @param penalties A list specifying the penalties to use. See details for
#'   more information.
#' @param simulation See documentation for AgeModel.
#' @param AgeVulnOffset Optional parameter for offset of age of
#'   vulnerability from age of maturity. A value of 0 means knife edge
#'   selectivity starts at 50% maturity (AgeMat), while -1 means 1 year
#'   before it. Thus (age of vulnerability = age of maturity +
#'   AgeVulnOffset). Defaults -1. Note that AgeMat is stochastic inside the
#'   function, so age at vulnerability will vary with it.
#' @param years Optional vector of years which is used by plotting
#'   functions. Defaults to 1:length(Catch).
#' @details The penalties list provides penalties and priors for up to four
#'   components of the analysis: initial depletion, carrying capacity (K),
#'   terminal year B/BMSY ('bstatus') and terminal year U/UMSY
#'   ('ustatus'). The user specifies which distribution to use and two
#'   distributional arguments for each of these metrics. For instance a
#'   ustatus~N(.5, .75) would be specified with list elements:
#'   ustatus.mean=0.5, ustatus.sd=0.75 and ustatus.dist=1. 'dist' options
#'   are 1=normal, 2=lognormal and 3=uniform. The uniform case passes
#'   arguments min and max, while the other two the mean and SD.
#'   deviation, and distribution type. Currently if no initial distribution
#'   is provided one is calculated internally but this needs to be
#'   revisiting and thus throws a warning.
#' @return A list containing depletion, SSB, and harvest rate (U) for
#'   posterior draws, and a vector of Keepers
#' @export
#'
run.SIR <- function(nrep, Catch, Taxon, penalties=NULL, Kprior=3,
                    Kscale=2, AgeVulnOffset=-1, years=NULL,
                    pct.keep=10, ProcessError=TRUE, simulation=NULL){
  check_penalties(penalties)
  pen <- penalties
  NY <- length(Catch)
  if(is.null(years)) years <- 1:NY
  stopifnot(length(Catch) == length(years))
  Bstore <- array(dim=c(NY,nrep))
  Dstore <- array(dim=c(NY,nrep))
  Ustore <- array(dim=c(NY,nrep))
  Bscaledstore <- Uscaledstore <- array(dim=c(NY,nrep))
  LikeStore <- rep(0, len=nrep)
  bmsy <- umsy <- cmsy <- crashed <- rep(NA, len=nrep)
  B <- array(dim=(NY+1))  #stock biomass
  ## Recruitment deviations (rows are replicates, columns years)
  recdevs <- matrix(NA, nrow=nrep, ncol=NY)
  ## Draw from biological  priors for SIR
  priors <- draw.priors(N=nrep, penalties=pen, Kscale=Kscale,
                        Kprior=Kprior, Catch=Catch, Taxon=Taxon)
  ## Add the intial conditions which are priors on K and depletion
  if(pen$initial.dist==2){
    ## lognormal
    InitialPrior <-
      rlnorm(nrep, meanlog=pen$initial.mean, sdlog=pen$initial.sd)
  } else {
    ## normal
    InitialPrior <-
      rnorm(nrep, mean=pen$initial.mean, sd=pen$initial.sd)
  }
  ## Carrying capacity (need to rename this at some point)
  if(is.null(pen$carry.dist)){
    message("No prior specified for K so generating one from max catch")
    ## catch info
    Cmax <- max(Catch,na.rm=TRUE)
    ## carrying capacity prior
    ## correlation between maximum catch and MSY
    Carry <- Kprior*Cmax
    Cprior <- runif(nrep,min=Carry/Kscale,max=Carry*Kscale)
    pen$carry.min <- Carry/Kscale; pen$carry.max <- Carry*Kscale;
    pen$carry.dist <- 3 # uniform
    check_penalties(pen)
  } else if (pen$carry.dist==1){
    Cprior <- rnorm(nrep, pen$carry.mean, pen$carry.sd)
  } else if(pen$carry.dist==2){
    Cprior <- rlnorm(nrep, pen$carry.mean, pen$carry.sd)
  }
  if(any(InitialPrior<0))
    warning('Some initial depletion values were negative')
  if(any(Cprior<0))
    warning('Some carrying capacity values were negative')
  ## Tack these onto the biological priors
  priors$draws$Cprior <- Cprior
  priors$draws$InitialPrior <- InitialPrior
  draws <- priors$draws
  ## Start of SIR
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
                    Weight,InitialDeplete,Sigma, AgeVulnOffset=AgeVulnOffset,
                    ProcessError=ProcessError,
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
      if(!is.null(pen$bstatus.dist)){
        x <- ifelse(pen$bstatus.dist==1, bscaled[NY], log(bscaled[NY]))
        loglike <- loglike +
          dnorm(x,mean=pen$bstatus.mean, sd=pen$bstatus.sd, log=TRUE)
      }
      ## If specified, include penalty for harvest rate in final year
      if(!is.null(pen$ustatus.dist)){
        x <- ifelse(pen$ustatus.dist==1, uscaled[NY], log(uscaled[NY]))
        loglike <- loglike +
          dnorm(x, mean=pen$ustatus.mean, sd=pen$ustatus.sd, log=TRUE)
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
  fit <- list(years=years, ssb=Bstore[,K], U=Ustore[,K], depletion=Dstore[,K],
              umsy=umsy[K], cmsy=cmsy[K], bmsy=bmsy[K], uscaled=Uscaledstore[,K],
              bscaled=Bscaledstore[,K], recdevs=recdevs, penalties=pen,
              Keepers=K, crashed=crashed, likes=LikeStore, Catch=Catch,
              draws=draws, Taxon=priors$Taxon)
  fit <- construct_srafit(fit)
  return(invisible(fit))
}
