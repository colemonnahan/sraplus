
#' Run the SIR algorithm for a draws and stock (temporary)
#' @param draws data.frame of prior draws for SIR
#' @return A list containing depletion, SSB, and harvest rate (U) for
#'   posterior draws, and a vector of Keepers
run.SIR <- function(draws, pct.keep=10){
  ## store results
  Bstore <- array(dim=c(NY,nrep))
  Dstore <- array(dim=c(NY,nrep))
  Ustore <- array(dim=c(NY,nrep))
  LikeStore <- array(dim=nrep)
  B <- array(dim=(NY+1))  #stock biomass
  ## #########################
  ## run iterations
  ## #########################
  nrep <- nrow(draws)
  for (irep in 1:nrep){   ##loop over replicates
    ## set priors
    InitialDeplete <- draws$InitialPrior[irep]
    Carry <- draws$Cprior[irep]
    ## life history
    Steep <- draws$h[irep]
    NatMort <- draws$M[irep]
    AgeMat <- round(draws$tm[irep],0)
    AgeMax <- ceiling(draws$tmax[irep])
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
    like <- dnorm(Deplete[NY],mean=FinalDepleteBest,sd=FinalDepleteCV)
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
              U=Ustore[,Keepers], Keepers=Keepers))
}
