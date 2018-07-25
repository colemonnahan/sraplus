### This file provides a quick demonstration of the package using a
### hard-coded example.
library(FishLife)
library(mvtnorm)
library(dplyr)
devtools::install('C:/Users/Cole/sraplus')
library(sraplus)
library(ggplot2)
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
## Run SIR to get posterior samples
fit <- run.SIR(nrep=nrep, Catch=Catch, Taxon=Taxon, InitialDepletePrior=.8,
               InitialDepleteCV=1, deplete.mean=2, deplete.cv=0.5,
               harvest.mean=0.5, harvest.sd=0.5,
               AgeVulnOffset=-2, years=2005:2015)

## Quick time series plots
par(mfrow=c(3,1))
plot_ssb(fit)
plot_bstatus(fit)
plot_ustatus(fit)

## plots of MSY reference points and NLL
plot_reference(fit)

## Look at biological prior vs posterior patterns
plot_draws(fit)

## Run a second fit and compare the differences
fit2 <- run.SIR(nrep=nrep, Catch=Catch, Taxon=Taxon, InitialDepletePrior=.8,
                InitialDepleteCV=1, deplete.mean=1.5, deplete.cv=0.3,
                harvest.mean=0.5, harvest.sd=0.1,
               AgeVulnOffset=-2, years=2005:2015)

## Compare two distinct fits (or more)
plot_terminal(fit, fit2)

## The return is an object of class 'srafit'
class(fit)
print(fit)
summary(fit)

## Plot management comparisons
axis.col <- gray(.5)
plot_fit(fit, fit2)
