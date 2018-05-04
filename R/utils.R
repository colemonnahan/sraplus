### Utility functions for package

#' Stochastic SRA to calculate the trajectory given parameters
#' @param Catch To do
AgeModel <- function(Catch,AgeMat, Steep,NatMort, AgeMax, Carry,Weight,InitialDeplete,Sigma) {
  ## process error devations
  devsvector=rnorm(NYears,mean=0,sd=Sigma)
  ## maximum age
  maxage=AgeMax
  ## natural survival
  surv = 1-NatMort
  ## 100% selectivity option
  if(AgeMat==1){vuln=rep(1,maxage)}
  ## assume knife edge selectivity
  if (AgeMat>1) {vuln=c(rep(0,AgeMat-1),rep(1,1+maxage-AgeMat))}
  ## assume logistic maturity
  mature = 1 / (1 + exp(AgeMat - 1:AgeMax))
  ## numbers at age
  num=rep(1,maxage)
  ## per recruit
  rzero=1  #to get sbpr
  for (a in 1:maxage) {
    ## equilibrium numbers at age
    if (a==1) num[a] = rzero   else num[a] = num[a-1]*surv
    ## plus group
    if (a==maxage) num[a] = num[a-1]*(surv/(1-surv))
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
  bio=c(1:NYears) ## spawning biomass
  pop=c(1:NYears) ## vulnerable biomass
  rec=c(1:NYears) ## recruitment
  hrstore=array(dim=NYears) ## harvest rate
  eggs = c(1:NYears)  ## spawning biomass
  tsurv=c(1:maxage) ## total survival
  for (y in 1:NYears) {
    ## vulnerable biomass
    Vpop=sum(num*vuln*Weight)
    ## harvest rate based on catch
    hr=min(.9,Catch[y]/Vpop)
    hrstore[y]=hr
    ## vulnerable biomass
    pop[y] = sum(num*vuln*Weight)
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
  } #end of loop over time
  return(pop)
} #end of function


#' Read catch data from file
#' @param File to do
ReadCatchData <- function(File){
  c="character"; n="numeric"
  D=read.csv(file=File,header=TRUE,colClasses=c(rep(c,6),rep(n,66)),sep=",")
  CatchData=D[,7:72]; colnames(CatchData)=seq(1950,2015)
  DD=list(StockInfo=D[,1:6],CatchData=CatchData)
  return(DD)
}


#' Read in stocks to consider and priors
#' @param File to do
ReadPrior <- function(File){
  c="character"; n="numeric"
  D=read.csv(file=File,header=TRUE,colClasses=c(rep(c,8),rep(n,6)),sep=",")
  return(D)
}
