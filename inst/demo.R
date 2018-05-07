## priors
## invertebrates
## figure out how far to push data sources on life history
## figure out priors for depletion - change in effort
	## know if species are caught by trawls?
## swept-area ratio - fraction between area trawled and continental shelf - correlates with fishing mortality rate
	## india - 10 - nothing long-lived is going to be around with that F


## identify what type of species from databases
  ## F prior
## how to use FMI as depletion measure
## relationship between trawl, effort, and F
## species x in area Z but not species y in area Z

### SHORT-TERM
## extend FishLife for invertebrates
## or use SeaLife/Base
## go from FMI to level of depletion
  ## FMI fishery score with more than 1 person
  ## related to level of depletion, at Bmsy 0.8 or higher, linear down to zero
## middle of April for invertebrates


rm(list=ls())

library(FishLife)
library(mvtnorm)
library(dplyr)
devtools::install('..')
library(sraplus)

###################
## directories
###################
main_dir <- file.path("C:\\merrill\\status_priors\\SSRA")
system.file('examples', 'simple', package='adnuts')

setwd(main_dir)

data_dir <- file.path(main_dir, "data")
R_dir <- file.path(main_dir, "R")

load("Return.Rdata")

###########################
## read global R functions
###########################

R_files <- list.files(R_dir)
ignore <- sapply(1:length(R_files), function(x) source(file.path(R_dir, R_files[x])))

###########################
## data
###########################
datafile <- file.path(data_dir, "Mexico FAO Catch Data.csv")
priorfile <- file.path(data_dir, "DeplPrior_wTaxa.csv")

## catch data
C <- ReadCatchData(File=datafile)
CatchData=C$CatchData
StockInfo=C$StockInfo
CStock=paste(StockInfo$Species..ASFIS.species.,StockInfo$Fishing.area..FAO.major.fishing.area.)

## prior info
Priors=ReadPrior(File=priorfile)
StockList=paste(Priors$Species..ASFIS.species.,Priors$Fishing.area..FAO.major.fishing.area.)

## indices
AllYears=seq(1950,2015)


###########################
## choose stock
###########################
iStock=3
jStock=which(CStock==StockList[iStock])

###########################
## FishLife
###########################

## settings
nrep=100000
Stock=StockList[iStock]

## initial depletion prior
InitialPrior=rnorm(nrep,Priors$InitialDepletePrior[iStock],Priors$InitialDepleteCV[iStock])

## catch info
Catch=as.numeric(CatchData[jStock,])
Cmax=max(Catch,na.rm=TRUE)

## final depletion prior
FinalDepleteBest=Priors$DepletePrior[iStock]
FinalDepleteCV=Priors$DepleteCV[iStock]

## carrying capacity prior
### correlation between maximum catch and MSY
Carry=Priors$Kprior[iStock]*Cmax
Cprior=runif(nrep,min=Carry/2,max=Carry*2)

#########################
## priors from FishLife
#########################
##------------ taxonomic information for stock -------------##
## draw taxonomic information from depletion prior file with stock info
Taxon <- Priors[iStock,which(colnames(Priors) %in% c("Class", "Order", "Family", "Genus", "Species"))]
Taxon[which(Taxon=="")] <- "predictive"

##------------ find information from FishLife -------------##
## find species within FishLife
sp <- Search_species(Class=Taxon[,"Class"], Order=Taxon[,"Order"], Family=Taxon[,"Family"], Genus=Taxon[,"Genus"], Species=Taxon[,"Species"], ParentChild_gz=Return$ParentChild_gz)$match_taxonomy

## find means and covariance matrix
Which <- grep(sp[1], Return$ParentChild_gz[,"ChildName"])
Mean <- Return$beta_gv[Which,]
Cov <- Return$Cov_gvv[Which,,]

##------------ deviates from multivariate normal -------------##
## parameters for multvariate normal distribution
params_mvn <- c("tm", "M", "tmax", "Loo", "K", "Winfinity") #"ln_var", "Loo", "h", "logitbound_h")
Mean_mvn <- Mean[which(names(Mean) %in% params_mvn)]
Cov_mvn <- Cov[which(rownames(Cov) %in% params_mvn), which(colnames(Cov) %in% params_mvn)]

## draw random deviates from multivariate normal distribution between
draws_mvn <- data.frame(rmvnorm(nrep, mean=Mean_mvn, sigma=Cov_mvn))

## exponentiate and add length-weight params
draws_mvn <- draws_mvn %>%
            mutate(Loo = exp(Loo)) %>%
            mutate(K = exp(K)) %>%
            mutate(Winfinity = exp(Winfinity)) %>%
            mutate(tmax = exp(tmax)) %>%
            mutate(tm = exp(tm)) %>%
            mutate(M = exp(M)) %>%
            mutate("lwb" = 3.04)  %>% ## from Froese, Thorson, and Reyes meta-analysis, mean value for all fish
            mutate("lwa" = Winfinity / (Loo ^ lwb))

##------------ recruitment deviations -------------##
Sigma <- rlnorm(nrep, Mean["ln_var"], 0.3)

##------------ recruitment deviations -------------##
Steep <- rlnorm(nrep, log(Mean["h"]), 0.1)


##------------ steepness -------------##
draws <- draws_mvn %>%
          mutate("Sigma" = Sigma) %>%
          mutate("h" = Steep)


## pairs plot of random draws
# pairs(draws)

## data inputs
Years=AllYears[which(is.na(Catch)==FALSE)]
Catch=Catch[which(is.na(Catch)==FALSE)]
NYears=length(Years)
NY=length(Years)

## store results
Bstore=array(dim=c(NY,nrep))
Dstore=array(dim=c(NY,nrep))
LikeStore=array(dim=nrep)
B=array(dim=(NY+1))  #stock biomass

###########################
## run iterations
###########################
for (irep in 1:nrep){   #loop over replicates
  # set priors
  InitialDeplete=InitialPrior[irep]
  Carry=Cprior[irep]

  ## life history
  Steep=draws$h[irep]
  NatMort=draws$M[irep]
  AgeMat=round(draws$tm[irep],0)
  AgeMax=ceiling(draws$tmax[irep])
  Sigma=draws$Sigma[irep]

  ## growth
  lwa <- draws$lwa[irep]
  lwb <- draws$lwb[irep]
  Linf <- draws$Loo[irep]
  vbk <- draws$K[irep]
  ages <- 1:AgeMax
  Length <- Linf*(1-exp(-vbk*ages))
  Weight <- lwa * Length ^ lwb

  pop=AgeModel(Catch, AgeMat, Steep,NatMort, AgeMax, Carry, Weight,InitialDeplete,Sigma)
  # plot(Years,pop,type="l",ylim=c(0,max(pop)))
  Deplete=pop/Carry

  like=dnorm(Deplete[NY],mean=FinalDepleteBest,sd=FinalDepleteCV)


  Bstore[,irep]=pop
  Dstore[,irep]=Deplete
  LikeStore[irep]=like
} #end of loop over replicates

## filter keepers
Nkeep=500
CumLike=c(0,cumsum(LikeStore))
BreakPoints=runif(Nkeep,min=0,max=CumLike[nrep+1])

Keepers=array(dim=(Nkeep))
k=1
for (i in 1:nrep){  #find the keepers
  Lower=CumLike[i]; Upper=CumLike[i+1]
  j=which(BreakPoints > Lower & BreakPoints < Upper)
  Nhits=length(j)
  if (Nhits>0) {Keepers[k:(k+Nhits-1)]= i; k=k+Nhits}

}

#plot the worm plot
plot(Years,Dstore[,Keepers[1]],ylim=c(0,2),type="l",xlab="Year",
     ylab="Depletion",main=Stock)

for (i in 2:Nkeep){
lines(Years,Dstore[,Keepers[i]],type="l")
}
lines(c(Years[1],Years[NY]),c(FinalDepleteBest,FinalDepleteBest),col="red",lwd=3)

