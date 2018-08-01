### This file provides a quick demonstration of the package using a
### hard-coded example.
library(FishLife)
library(mvtnorm)
library(dplyr)
## devtools::document('..')
## devtools::install('C:/Users/Cole/sraplus', quick=TRUE)
devtools::load_all('C:/Users/Cole/sraplus')
library(sraplus)
library(ggplot2)
## You need to put the file "Return.RData" in the "data" folder of the
## package. It's too large to include it and should be removed later. Comes
## from FishLife
data(Return)

## A simulated stock history for demonstration purposes
nrep <- 10000 # total reps, 10% will be kept
set.seed(2323)
Catch <- runif(11, 50000, 300000)
Taxon <- c(Class="Actinopterygii", Order="Perciformes",
           Family="Scombridae", Genus="Thunnus", Species="albacares")
## Define the penalties. Defaults to normal distribution and a initial
## distribution that is based on max catch (needs to be updated)
pen <- list(bstatus.mean=0, bstatus.sd=0.5, bstatus.dist=2,
            ustatus.mean=0, ustatus.sd=0.25, ustatus.dist=2,
            initial.mean=0, initial.sd=.3, initial.dist=2)
## Run SIR to get posterior samples
fit <- run.SIR(nrep=nrep, Catch=Catch, Taxon=Taxon, penalties=pen,
               AgeVulnOffset=-2, years=2005:2015)

## Quick time series plots
par(mfrow=c(3,1))
plot_ssb(fit)
plot_bstatus(fit)
plot_ustatus(fit)

## Look at biological prior vs posterior patterns, this documentation and
## function need updating. Red points are crashed, black points not kept,
## and green kept.
plot_draws(fit)

## Run an arbitrary second fit and compare the differences.
pen$ustatus.mean <- -.5
pen$bstatus.sd <- .25
fit2 <- run.SIR(nrep=nrep, Catch=Catch, Taxon=Taxon,
                penalties=pen, AgeVulnOffset=-2, years=2005:2015)

## Compare two distinct fits (or more)
plot_terminal(fit, fit2)
## plots of MSY reference points and NLL
plot_reference(fit, fit2)

## The return is an object of class 'srafit'. Not heavily used for now but
## probably will be useful later.
class(fit)
print(fit)
summary(fit) # the same for now

## Plot management comparisons using specially-designed function.
plot_fit(fit, fit2)

