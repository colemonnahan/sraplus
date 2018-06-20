#' Stochastic SRA to calculate the trajectory given parameters
#' @param Catch To do
#' @param AgeMat Age at maturity
#' @param Steep Steepness
#' @param NatMort Natural mortality
#' @param AgeMax Maximum age, above which there is a plus group.
#' @param Carry Carrying capacity.
#' @param Weight Vector of weights at age, calculated from growth
#'   parameters a and b.
#' @param InitialDeplete Initial level of depletion.
#' @param Sigma Process error variance, i.e. N(0,Sigma).
#' @param ProcessError A flag for whether to include process error,
#'   defaults to TRUE.
#' @export
AgeModel <- function(Catch, AgeMat, Steep, NatMort, AgeMax,
                     Carry, Weight, InitialDeplete, Sigma,
                     ProcessError=TRUE, simulation=NULL){
  stopifnot(AgeMat>0)
  use.sim <- !is.null(simulation)
  if(use.sim){
    stopifnot(is.numeric(simulation$FishMort))
    stopifnot(simulation$FishMort>=0); stopifnot(simulation$FishMort<=1)
    NYears <- simulation$NYears
  } else {
    NYears <- length(Catch)
  }
  ## The process error deviations
  if(ProcessError){
    devs <- rnorm(NYears,mean=0,sd=Sigma)
  } else {
    devs <- rep(0, len=NYears)
    Sigma <- 0 # since used below in correction
  }
  ## maximum age
  maxage <- AgeMax
  ## natural survival
  surv  <-  1-NatMort
  ## num is numbers at age; vuln is vulnerability
  num <- vuln <- rep(1,maxage)
  ## 100% selectivity option
  if(AgeMat>1){
    ## assume knife edge selectivity
    vuln[1:(AgeMat-1)] <- 0
  }
  ## assume logistic maturity with 50% at AgeMat
  mature  <-  1 / (1 + exp(AgeMat - 1:AgeMax))
  ## per recruit
  rzero <- 1  #to get sbpr
  for (a in 1:maxage) {
    ## equilibrium numbers at age
    if (a==1) num[a]  <-  rzero   else num[a]  <-  num[a-1]*surv
    ## plus group
    if (a==maxage) num[a]  <-  num[a-1]*(surv/(1-surv))
  }
  ## equilibrium spawning biomass
  bzeroInit <- sum(Weight*num*mature)# biomass per recruit
  ## rzero scaled to carrying capacity
  rzero <- Carry/bzeroInit
  bzero <- bzeroInit*rzero
  ## include process error
  rzeroC  <-  rzero*exp(-Sigma*2/2) #correct for recruitment devs
  ## Beverton-Holt stock-recruit parameters
  alpha <- (bzero*(1-Steep))/(4*Steep*rzeroC)
  betav <- (5*Steep-1)/(4*Steep*rzeroC)
  ## set initial numbers at age
  num0 <- num <- num*rzero*InitialDeplete
  ##  now loop over time
  ## bio is spawning biomass; pop is vulnerable biomass; rec is recruits;
  ## hrstore is the harvest rate
  hrstore <-  bio <- pop <- rec <- eggs <-
    rep(NA, length=NYears)
  ##  tsurv <- c(1:maxage) ## total survival
  for (y in 1:NYears) {
    ## vulnerable biomass
    Vpop <- sum(num*vuln*Weight)
    ## if simulating a population, calculate Catch
    if(use.sim)
      Catch[y] <- simulation$FishMort*Vpop
    ## harvest rate based on catch
    hr <- Catch[y]/Vpop
    ## hr <- min(.9,Catch[y]/Vpop)
    hrstore[y] <- hr
    if(hr>1) break ## break if caught more than available
    ## vulnerable biomass
    pop[y]  <-  sum(num*vuln*Weight)
    if(pop[y]<0) break ## break if negative
    ## spawning biomass
    eggs[y]  <-  sum(num*mature*Weight)
    ## recruits with process error
    rec[y]  <-  eggs[y]/(alpha+betav*eggs[y])  * exp(devs[y])
    ## abundance in plus group
    num[maxage] <- num[maxage]*surv*(1-hr*vuln[maxage])+num[maxage-1]*surv*(1-hr*vuln[maxage-1])  #update plus group
    ## total survival
    tsurv <- (1-vuln*hr)*surv                   #vector of total survivals
    ## fill in age classes
    num[2:(maxage-1)] <- num[1:(maxage-2)]*tsurv[1:(maxage-2)]        #update age classes 2 to maxage-1
    if (y==1) {
      num[1] <- num[1]
    } else {
      ## recruitment
      num[1] <-  rec[y-1]
    }
    if(any(num<0)) {message("breaking due to crashed pop") ;break}
  } #end of loop over time

  ## Calculate UMSY
  ##  Weight <- c(0.35, 0.57, 0.78, 0.97, 1.14, 1.27, 1.37, 1.46, 1.52, 1.57, 1.60, 1.63, 1.65, 1.67, 1.68, 1.68)
  get.equilibrium.catch <- function(U){
    ## equilibrium numbers at age under rate U
    num[1]  <-  rzeroC
    for (a in 2:(maxage-1))
      num[a]  <-  (1-vuln[a-1]*U)*num[a-1]*surv
    ## plus group
    num[maxage]  <- num[a-1]*( (1-vuln[a-1]*U)*surv / (1-(1-vuln[a-1]*U)*surv) )
    ## Now get equilibrium recruitment for U
    ## SBPR0 <- sum(Weight*mature*num0)/rzeroC
    ## SSB0 <- SBPR0*rzeroC
    ## spwaning biomass per recruit fishing at U
    SBPR <- sum(Weight*mature*num)/rzeroC
    ## recruits in equilibrium fishing at U
    R <- max(0, (SBPR-alpha)/(betav*SBPR))
    ## and equlibrium catch per recruit?
    Ca <- sum(num*vuln*Weight*U)
    YPR <- Ca/rzeroC
    ## Calculate catch for all recruits
    return(R*YPR)
  }
  fit <- optimize(get.equilibrium.catch, interval=c(0,1), maximum=TRUE)
  out <- list(pop=pop, hr=hrstore, umsy=fit$maximum, cmsy=fit$objective)
  if(use.sim){
    u.seq <- seq(0,1, len=100)
    c.seq <- sapply(u.seq, function(u) get.equilibrium.catch(u))
    o <- list(num0=num0, num=num, Catch=Catch,
              CatchEquilibrium=get.equilibrium.catch(simulation$FishMort),
              u.seq=u.seq, c.seq=c.seq)
    out <- c(out, o)
  }
  return(out)
} #end of function

