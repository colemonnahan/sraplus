#' Stochastic SRA to calculate the trajectory given parameters
#' @param Catch To do
#' @export
AgeModel <- function(Catch,AgeMat, Steep,NatMort, AgeMax, Carry,Weight,InitialDeplete,Sigma) {
  ## process error devations
  devsvector <- rnorm(NYears,mean=0,sd=Sigma)
  ## maximum age
  maxage <- AgeMax
  ## natural survival
  surv  <-  1-NatMort
  ## 100% selectivity option
  if(AgeMat==1){vuln <- rep(1,maxage)}
  ## assume knife edge selectivity
  if (AgeMat>1) {vuln <- c(rep(0,AgeMat-1),rep(1,1+maxage-AgeMat))}
  stopifnot(AgeMat>0)
  ## assume logistic maturity
  mature  <-  1 / (1 + exp(AgeMat - 1:AgeMax))
  ## numbers at age
  num <- rep(1,maxage)
  ## per recruit
  rzero <- 1  #to get sbpr
  for (a in 1:maxage) {
    ## equilibrium numbers at age
    if (a==1) num[a]  <-  rzero   else num[a]  <-  num[a-1]*surv
    ## plus group
    if (a==maxage) num[a]  <-  num[a-1]*(surv/(1-surv))
  }
  ## equilibrium spawning biomass
  bzeroInit=sum(Weight*num*mature)# biomass per recruit
                                        # rzero scaled to carrying capacity
  rzero=Carry/bzeroInit
  bzero=bzeroInit*rzero
  ## include process error
  rzeroC = rzero*exp(-Sigma*2/2) #correct for recruitment devs
  ## Beverton-Holt stock-recruit parameters
  alpha=(bzero*(1-Steep))/(4*Steep*rzeroC)
  betav=(5*Steep-1)/(4*Steep*rzeroC)
                                        # set initial numbers at age
  num=num*rzero*InitialDeplete
  ##  now loop over time
  bio=rep(NA, length=NYears) ## spawning biomass
  pop=rep(NA, length=NYears) ## vulnerable biomass
  rec=rep(NA, length=NYears) ## recruitment
  hrstore=array(dim=NYears) ## harvest rate
  eggs = rep(NA, length=NYears)  ## spawning biomass
  tsurv=c(1:maxage) ## total survival
  for (y in 1:NYears) {
    ## vulnerable biomass
    Vpop=sum(num*vuln*Weight)
    ## harvest rate based on catch
    hr <- Catch[y]/Vpop
    ## hr=min(.9,Catch[y]/Vpop)
    hrstore[y]=hr
    if(hr>1) break ## break if caught more than available
    ## vulnerable biomass
    pop[y] = sum(num*vuln*Weight)
    if(pop[y]<0) break ## break if negative
    ## spawning biomass
    eggs[y] = sum(num*mature*Weight)
    ## process error
    dev=devsvector[y]
    ## recruits with process error
    rec[y] = eggs[y]/(alpha+betav*eggs[y])  * exp(dev)
    ## abundance in plus group
    num[maxage]=num[maxage]*surv*(1-hr*vuln[maxage])+num[maxage-1]*surv*(1-hr*vuln[maxage-1])  #update plus group
    ## total survival
    tsurv=(1-vuln*hr)*surv                   #vector of total survivals
    ## fill in age classes
    num[2:(maxage-1)]=num[1:(maxage-2)]*tsurv[1:(maxage-2)]        #update age classes 2 to maxage-1
                                        # print(num[1:10])
    if (y==1) {num[1]=num[1]} else {num[1]= rec[y-1] }  #recruitment
    if(any(num<0)) break
  } #end of loop over time
  return(list(pop=pop, hr=hrstore))
} #end of function

