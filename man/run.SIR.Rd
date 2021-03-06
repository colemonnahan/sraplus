% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SIR.R
\name{run.SIR}
\alias{run.SIR}
\title{Run the SIR algorithm for a draws and stock (temporary)}
\usage{
run.SIR(nrep, Catch, Taxon, penalties = NULL, Kprior = 3, Kscale = 2,
  AgeVulnOffset = -1, years = NULL, pct.keep = 10, ProcessError = TRUE,
  simulation = NULL)
}
\arguments{
\item{nrep}{Number of replicates to run.}

\item{Catch}{Vector of catches, one for each year.}

\item{Taxon}{A list with taxonomic information, including 'Class',
'Order', 'Family', 'Genus', 'Species'. Any blanks will be converted to
'predictive' internally when querying FishLife.}

\item{penalties}{A list specifying the penalties to use. See details for
more information.}

\item{Kprior}{Maybe a multiplier on K?}

\item{Kscale}{Scalar to control the initial biomass which is generated}

\item{AgeVulnOffset}{Optional parameter for offset of age of
vulnerability from age of maturity. A value of 0 means knife edge
selectivity starts at 50% maturity (AgeMat), while -1 means 1 year
before it. Thus (age of vulnerability = age of maturity +
AgeVulnOffset). Defaults -1. Note that AgeMat is stochastic inside the
function, so age at vulnerability will vary with it.}

\item{years}{Optional vector of years which is used by plotting
functions. Defaults to 1:length(Catch).}

\item{pct.keep}{The percentage of "keepers" from the total. Default
10\%.}

\item{ProcessError}{Flag for whether to include process error
deviations, passed on to model. Defaults to TRUE.}

\item{simulation}{See documentation for AgeModel.}
}
\value{
A list containing depletion, SSB, and harvest rate (U) for
  posterior draws, and a vector of Keepers
}
\description{
Run the SIR algorithm for a draws and stock (temporary)
}
\details{
The penalties list provides penalties and priors for up to four
  components of the analysis: initial depletion, carrying capacity (K),
  terminal year B/BMSY ('bstatus') and terminal year U/UMSY
  ('ustatus'). The user specifies which distribution to use and two
  distributional arguments for each of these metrics. For instance a
  ustatus~N(.5, .75) would be specified with list elements:
  ustatus.mean=0.5, ustatus.sd=0.75 and ustatus.dist=1. 'dist' options
  are 1=normal, 2=lognormal and 3=uniform. The uniform case passes
  arguments min and max, while the other two the mean and SD.
  deviation, and distribution type. Currently if no initial distribution
  is provided one is calculated internally but this needs to be
  revisiting and thus throws a warning.
}
