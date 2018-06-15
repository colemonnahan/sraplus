#' Stochastic SRA to calculate the trajectory given parameters
#' @param Catch To do
#' @export
AgeModel <- function(Catch, AgeMat, Steep, NatMort, AgeMax,
                     Carry, Weight, InitialDeplete, Sigma) {
  ## process error devations
  stopifnot(AgeMat>0)
  NYears <- length(Catch)
  devsvector <- rnorm(NYears,mean=0,sd=Sigma)
  ## maximum age
  maxage <- AgeMax
  ## natural survival
  surv  <-  1-NatMort
  ## num is numbers at age; vuln is vulnerability
  num <- vuln <- rep(1,maxage)
  ## 100% selectivity option
  if(AgeMat>1){
    ## assume knife edge selectivity
    vuln[1:(AgeMat-1)] <- 1
  }
  ## assume logistic maturity
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
  num <- num*rzero*InitialDeplete
  ##  now loop over time
  ## bio is spawning biomass; pop is vulnerable biomass; rec is recruits;
  ## hrstore is the harvest rate
  hrstore <-  bio <- pop <- rec <- eggs <-
    rep(NA, length=NYears)
  ##  tsurv <- c(1:maxage) ## total survival
  for (y in 1:NYears) {
    ## vulnerable biomass
    Vpop <- sum(num*vuln*Weight)
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
    rec[y]  <-  eggs[y]/(alpha+betav*eggs[y])  * exp(devsvector[y])
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
  return(list(pop=pop, hr=hrstore))
} #end of function

