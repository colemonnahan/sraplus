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
devtools::load_all('C:/Users/Cole/sraplus')
library(sraplus)

###################
## directories
###################
main_dir <- system.file(package='sraplus')
setwd(main_dir)
data_dir <- file.path(main_dir, "examples")

###########################
## data
###########################
load(file.path(data_dir, "Return.RData"))
datafile <- file.path(data_dir, "Mexico FAO Catch Data.csv")
priorfile <- file.path(data_dir, "DeplPrior_wTaxa.csv")

## catch data
C <- ReadCatchData(File=datafile)
CatchData <- C$CatchData
StockInfo <- C$StockInfo
CStock <- paste(StockInfo$Species..ASFIS.species.,StockInfo$Fishing.area..FAO.major.fishing.area.)

## prior info
Priors <- ReadPrior(File=priorfile)
StockList <- paste(Priors$Species..ASFIS.species.,Priors$Fishing.area..FAO.major.fishing.area.)

## indices
AllYears=seq(1950,2015)


###########################
## choose stock
###########################
iStock <- 1
jStock <- which(CStock==StockList[iStock])
## data inputs for chosen stock
Catch <- as.numeric(CatchData[jStock,])
Years <- AllYears[which(is.na(Catch)==FALSE)]
Catch <- Catch[which(is.na(Catch)==FALSE)]
NYears <- length(Years)
NY <- length(Years)
## final depletion prior
FinalDepleteBest <- Priors$DepletePrior[iStock]
FinalDepleteCV <- Priors$DepleteCV[iStock]

## settings
nrep <- 1000
set.seed(2323)
## Draw from priors for SIR
draws <- draw.priors(nrep, iStock)
## Run SIR to get posterior samples
temp <- run.SIR(draws)
tail(temp$depletion[1,])
## [1] 0.9965566 1.0480772 1.5582368 0.7907136 1.4609847 1.1304281

## Quick exploratory plots
par(mfrow=c(1,2))
## Plot depletion trajectories
plot(Years, Years, ylim=c(0,2),type="n",xlab="Year",
     ylab="Depletion",main=StockList[iStock])
trash <- apply(temp$depletion, 2, function(i) lines(Years,y=i))
abline(h=FinalDepleteBest,col="red",lwd=3)
## Same but for biomass
plot(Years,Years, ylim=c(0,1.1*max(temp$ssb)), type="n",xlab="Year",
     ylab="Vulnerable Biomass",main=StockList[iStock])
trash <- apply(temp$ssb, 2, function(i) lines(Years,y=i))

## Pairs plot of parameters priors and posterior
dev.new()
post.ind <- 1:ncol(temp$depletion) %in% temp$Keepers # which prior draws were "kept"
pairs(draws[,1:5], col=ifelse(post.ind, 'red', 'black'),
      cex=ifelse(post.ind, 1, .1))
