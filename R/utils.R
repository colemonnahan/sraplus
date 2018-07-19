### Utility functions for package

#' Plot spawning biomass time series for posterior draws with catches.
#' @template plot_args
#' @param ylim The ylim to use, if not specified will be calculated
#'   internally.
#' @export
plot_ssb <- function(fit, ylim=NULL){
  if(is.null(fit$year)) fit$year <- 1:nrow(fit$ssb)
  if(is.null(ylim)) ylim <- c(0, 1.05*max(fit$ssb))
  plot(x=fit$year, y=fit$year, ylim=ylim, type="n",xlab='Year',
       ylab="Vulnerable Biomass", main=NA)
  ## Add individual trajectories that were kept
  trash <- apply(fit$ssb, 2, function(i)
    lines(fit$year, y=i, col=rgb(0,0,0,.1)))
  ## Add catch
  points(fit$year, fit$Catch, type='h', col='blue', lwd=3)
}

#' Plot the prior draws from a fit, showing which crashed the population
#' and which were kept
#' @template plot_args
#' @export
plot_draws <- function(fit){
  n <- nrow(fit$draws)
  col <- rep('black', len=n)
  col[fit$crashed] <- 'red'
  col[unique(fit$Keepers)] <- 'green'
  ## want the order to be black red green I think
  col <- factor(col, levels=c('red', 'black', 'green'))
  ind <- order(col)
  pairs(fit$draws[ind,], col=col[ind], upper.panel=NULL, pch='.')
}

#' Plot the realized recruitment deviations from a fit, showing which
#' crashed the population and which were kept
#' @template plot_args
#' @export
plot_recdevs <- function(fit){
  n <- nrow(fit$draws)
  col <- rep('black', len=n)
  col[fit$crashed] <- 'red'
  col[unique(fit$Keepers)] <- 'green'
  ## want the order to be black red green I think
  col <- factor(col, levels=c('red', 'black', 'green'))
  ind <- order(col)
  col <- col[ind]
  recdevs <- fit$recdevs[ind,]
  years <- 1:ncol(recdevs)
  ylim <- c(-1,1)*1.05*max(abs(recdevs))
  plot(years,  ylim=ylim, type="n",xlab=NA, pch='.',
       ylab="Recruitment deviation")
  for(i in 1:nrow(recdevs)){
    points(jitter(years), y=recdevs[i,], col=as.character(col[i]))
  }
}

