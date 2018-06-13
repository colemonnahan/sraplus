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
iStock <- 4
jStock <- which(CStock==StockList[iStock])
# data inputs for chosen stock
Catch <- as.numeric(CatchData[jStock,])
Years <- AllYears[which(is.na(Catch)==FALSE)]
Catch <- Catch[which(is.na(Catch)==FALSE)]
NYears <- length(Years)
NY <- length(Years)

## final depletion prior
CStock <- "Marine fishes nei Pacific, Eastern Central"
Years <- 1956:2015
FinalDepleteBest <- Priors$DepletePrior[iStock]
FinalDepleteCV <- Priors$DepleteCV[iStock]
Kprior <- Priors$Kprior[iStock]
InitialDepletePrior <- Priors$InitialDepletePrior[iStock]
InitialDepleteCV <- Priors$InitialDepleteCV[iStock]
Taxon <- c(Class="Actinopterygii", Order="", Family="", Genus="", Species="")
  Priors[iStock, c("Class", "Order", "Family", "Genus", "Species")]
Years <- 1950:2015
NYears <- length(Years)
Catch <- c(11200, 11300, 8700, 10000, 5800, 3400, 5200, 2700, 19200, 16400,
           25600, 40000, 31800, 34000, 30000, 38400, 50000, 78000, 88100, 86500,
           91500, 111100, 117100, 122300, 18936, 25178, 40793, 2575, 97035, 157754,
           137998, 247485, 32678, 43009, 93478, 107547, 111452, 90206, 85661, 78127,
           243749, 269734, 147620, 120349, 174176, 238947, 180690, 131611, 113369,
           66776, 83729, 71432, 74800, 72603, 59766, 72231, 75397, 51974, 59222,
           46688, 48924, 47367, 78814, 18464, 22498, 32557)

## settings
devtools::load_all('C:/Users/Cole/sraplus')
library(sraplus)
nrep <- 10000
set.seed(2323)
## Draw from priors for SIR
draws <- draw.priors(N=nrep, InitialDepletePrior, InitialDepleteCV,
                     Kprior, Catch, Taxon)
## Run SIR to get posterior samples
fit1 <- run.SIR(NY=NY, Catch=Catch, draws=draws, deplete.mean=FinalDepleteBest,
                deplete.cv=FinalDepleteCV )

