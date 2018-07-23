### This file provides a quick demonstration of the package using a
### hard-coded example.
library(FishLife)
library(mvtnorm)
library(dplyr)
devtools::install('C:/Users/Cole/sraplus')
library(sraplus)
## You need to put the file "Return.RData" in the "data" folder of the
## package. It's too large to include it and should be removed later. Comes
## from FishLife
data(Return)

## A simulated stock history for demonstration purposes

## settings
nrep <- 10000
set.seed(2323)
Catch <- runif(11, 50000, 300000)
Taxon <- c(Class="Actinopterygii", Order="Perciformes",
           Family="Scombridae", Genus="Thunnus", Species="albacares")
## Draw from priors for SIR
draws <- draw.priors(N=nrep, InitialDepletePrior=0.8,
                     InitialDepleteCV=0.2,
                     Kprior=3, Catch=Catch, Taxon=Taxon)
## Run SIR to get posterior samples
fit <- run.SIR(Catch=Catch, draws=draws, deplete.mean=2,
                deplete.cv=0.5)

## Quick plot to show the output
par(mfrow=c(3,1))
fit$years <- 2005:2015
plot_ssb(fit)
plot_bstatus(fit)
plot_ustatus(fit)
## plots of MSY reference points and NLL
plot_reference(fit)
