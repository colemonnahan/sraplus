% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot_draws}
\alias{plot_draws}
\title{Plot the prior draws from a fit, showing which crashed the population
and which were kept}
\usage{
plot_draws(fit)
}
\arguments{
\item{fit}{A returned list from run.SIR}

\item{...}{A series of fits as returned by run.SIR}

\item{names}{An optional character string with length equal to the
number of fits, used to identify fits on the plot.}

\item{plot}{Boolean for whether to print the plot.}
}
\value{
Functions using base plots return nothing, but those using
  ggplot2 return an invisible ggplot object which can then be used to
  print or save as normally done. In all cases if 'plot=TRUE' then a
  plot is created.
}
\description{
Plot the prior draws from a fit, showing which crashed the population
and which were kept
}
