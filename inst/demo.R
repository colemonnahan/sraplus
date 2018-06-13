### This file provides a quick demonstration of the package using a
### hard-coded example.
library(FishLife)
library(mvtnorm)
library(dplyr)
devtools::load_all('C:/Users/Cole/sraplus')
library(sraplus)
## You need to put the file "Return.RData" in the "data" folder of the
## package. It's too large to include it and should be removed later. Comes
## from FishLife
data(Return)
str(Return)

## Hard-coded values for a specific stock
CStock <- "Marine fishes nei Pacific, Eastern Central"
Years <- 1956:2015
FinalDepleteBest <- 0.15
FinalDepleteCV <- 0.2
Kprior <- 5
InitialDepletePrior <- 0.8
InitialDepleteCV <- 0.2
Taxon <- c(Class="Actinopterygii", Order="", Family="", Genus="",
           Species="")
Years <- 1950:2015
NY <- length(Years)
Catch <- c(11200, 11300, 8700, 10000, 5800, 3400, 5200, 2700, 19200, 16400,
           25600, 40000, 31800, 34000, 30000, 38400, 50000, 78000, 88100, 86500,
           91500, 111100, 117100, 122300, 18936, 25178, 40793, 2575, 97035, 157754,
           137998, 247485, 32678, 43009, 93478, 107547, 111452, 90206, 85661, 78127,
           243749, 269734, 147620, 120349, 174176, 238947, 180690, 131611, 113369,
           66776, 83729, 71432, 74800, 72603, 59766, 72231, 75397, 51974, 59222,
           46688, 48924, 47367, 78814, 18464, 22498, 32557)

## settings
nrep <- 10000
set.seed(2323)
## Draw from priors for SIR
draws <- draw.priors(N=nrep, InitialDepletePrior, InitialDepleteCV,
                     Kprior, Catch, Taxon)
## Run SIR to get posterior samples
fit <- run.SIR(Catch=Catch, draws=draws, deplete.mean=FinalDepleteBest,
                deplete.cv=FinalDepleteCV )

## Quick plot to show the output
par(mfrow=c(3,1))
plot(Years, Years, ylim=c(0,2),type="n",xlab=NA,
     ylab="Depletion",main=CStock)
trash <- apply(fit$depletion, 2, function(i) lines(Years,y=i, col=rgb(0,0,0,.1)))
plot(Years,Years, ylim= c(0,1.05*max(fit$ssb)), type="n",xlab=NA,
     ylab="Vulnerable Biomass")
trash <- apply(fit$ssb, 2, function(i) lines(Years,y=i, col=rgb(0,0,0,.1)))
## Same but for harvest rate (U)
plot(Years,Years, ylim=c(0,1), type="n",xlab="Year",
     ylab="Harvest rate (U)")
trash <- apply(fit$U, 2, function(i) lines(Years,y=i, col=rgb(0,0,0,.1)))
