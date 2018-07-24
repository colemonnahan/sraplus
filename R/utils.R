### Utility functions for package

#' Plot posterior histograms of MSY reference points and the log posterior
#' density
#'
#' @template plot_args
#' @export
plot_reference <- function(fit){
  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))
  par(mfrow=c(2,2), mgp=c(1.1, .3, 0), tck=-.02, mar=c(2.5,2.5,.5,.5))
  hist(log(fit$likes[fit$likes>0]), main=NA, xlab="Likelihood", ylab=NA)
  hist(fit$bmsy, main=NA, xlab="BMSY", ylab=NA)
  hist(fit$umsy, main=NA, xlab="UMSY", ylab=NA, xlim=c(0,1))
  hist(fit$cmsy, main=NA, xlab="MSY", ylab=NA)
}


#' Plot B/BMSY (biomass status) time series for posterior draws, with
#'   terminal year penalty shown if used.
#' @template plot_args
#' @param ylim The ylim to use, if not specified will be calculated
#'   internally.
#' @export
plot_bstatus <- function(fit, ylim=NULL){
  if(is.null(fit$year)) fit$year <- 1:nrow(fit$bscaled)
  if(is.null(ylim)) ylim <- c(0, 1.05*max(fit$bscaled))
  plot(x=fit$year, y=fit$year, ylim=ylim, type="n",xlab='Year',
       ylab="B/BMSY")
  trash <- apply(fit$bscaled, 2, function(i)
    lines(x=fit$year, y=i, col=rgb(0,0,0,.1)))
  ci <- qnorm(c(0.025, .975), mean=fit$penalties$deplete.mean, fit$penalties$deplete.cv)
  lines(x=rep(tail(fit$year,1),2), ci, col=2, lwd=2)
  points(x=tail(fit$year,1), y=fit$penalties$deplete.mean, col=2, cex=1.5,
         pch=16)
  abline(h=1, col='red', lty=3, lwd=2)
}

#' Plot U/UMSY (exploitation status) time series for posterior draws, with
#'   terminal year penalty shown if used.
#' @template plot_args
#' @param ylim The ylim to use, if not specified will be calculated
#'   internally.
#' @export
plot_ustatus <- function(fit, ylim=NULL){
  if(is.null(fit$year)) fit$year <- 1:nrow(fit$uscaled)
  if(is.null(ylim)) ylim <- c(0, 1.05*max(fit$uscaled))
  plot(x=fit$year, y=fit$year, ylim=ylim, type="n",xlab='Year',
       ylab="U/UMSY")
  trash <- apply(fit$uscaled, 2, function(i) lines(fit$year,y=i, col=rgb(0,0,0,.1)))
  ci <- qnorm(c(0.025, .975), mean=fit$penalties$deplete.mean, fit$penalties$deplete.cv)
  lines(x=rep(tail(fit$year,1),2), ci, col=2, lwd=2)
  points(x=tail(fit$year,1), y=fit$penalties$deplete.mean, col=2, cex=1.5,
         pch=16)
  abline(h=1, col='red', lty=3, lwd=2)
}


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

