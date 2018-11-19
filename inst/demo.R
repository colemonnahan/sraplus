  ### This file provides a quick demonstration of the package using a
  ### hard-coded example. In development and subject to change!

  ## devtools::document('..')
  ## devtools::install('C:/Users/Cole/sraplus', quick=TRUE)
  ## devtools::load_all('C:/Users/Cole/sraplus')
  # devtools::install_github(repo='colemonnahan/sraplus', quick=TRUE)
  library(FishLife)
  library(mvtnorm)
  library(dplyr)
  library(sraplus)
  library(ggplot2)
  ## You need to put the file "Return.RData" in the "data" folder of the
  ## package. It's too large to include it and should be removed later. Comes
  ## from FishLife
  data(Return)

  ## A simulated stock history for demonstration purposes
  nrep <- 20000 # total reps, 10% will be kept
  set.seed(2323)
  Catch <- runif(11, 50000, 300000)
  Taxon <- c(Class="Actinopterygii", Order="Perciformes",
             Family="Scombridae", Genus="Thunnus", Species="albacares")
  ## Define the penalties. Defaults to a uniform carrying capacity
  ## distribution that is based on max catch (needs to be updated). bstatus =
  ## terminal B/BMSY; ustatus=terminal U/UMSY; initial=initial depletion
  pen <- list(bstatus.mean=0, bstatus.sd=0.5, bstatus.dist=2,
              ustatus.mean=0.5, ustatus.sd=0.25, ustatus.dist=2,
              initial.mean=0, initial.sd=.3, initial.dist=2)
  ## Run SIR to get posterior samples
  fit <- run.SIR(nrep=nrep, Catch=Catch, Taxon=Taxon, penalties=pen,
                  years=2005:2015, inst_f = FALSE)

## Quick time series plots
par(mfrow=c(3,1))
plot_ssb(fit)
plot_bstatus(fit)
plot_ustatus(fit)

## Look at biological prior vs posterior patterns, this documentation and
## function need updating. Red points are crashed, black points not kept,
## and green kept.
plot_draws(fit)
plot_fit(fit)
## Run an arbitrary second fit and compare the differences.



pen2 <- list(bstatus.mean=log(1), bstatus.sd=0.5, bstatus.dist=2,
            ustatus.mean= log(1.5), ustatus.sd=0.25, ustatus.dist=2,
            initial.mean=log(0.5), initial.sd=.3, initial.dist=2 ,
            ## now we specify a lognormal distribution for K explicitly
            carry.mean=13.8, carry.sd=.2, carry.dist=2)
## Also change the AgeVulnOffset from default of -1
fit2 <- run.SIR(nrep=nrep, Catch=Catch, Taxon=Taxon,
                penalties=pen2, years=2005:2015, AgeVulnOffset = -1)

## Compare prior and posterior for two distinct fits (or more)
plot_penalties(fit, fit2)
## plots of MSY reference points and NLL
plot_reference(fit, fit2)

## The return is an object of class 'srafit'. Not heavily used for now but
## probably will be useful later.
class(fit)
print(fit)
## credible intervals for prior and posterior on key management targets
summary(fit)

## Plot management comparisons using specially-designed function.
plot_fit(fit, fit2)

