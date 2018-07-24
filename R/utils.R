### Utility functions for package

#' Plot terminal year U/UMSY and B/BMSY posteriors and priors (if
#' applicable) for a series of fits
#' @param ... A series of fits as returned by run.SIR
#' @param names An optional character string with length equal to the
#' number of fits, used to identify fits on the plot.
#' @param plot Boolean for whether to print the plot.
#' @param xlim Optional xlim to override the defaults chosen by ggplot
#' @param q Optional quantiles between which to calculate prior
#' density. Defaults to c(0.005,0.995), that is 99\%.
#' @return An invisible ggplot2 object which can be printed or saved.
#' @export
plot_terminal <- function(..., names=NULL, plot=TRUE, xlim=NULL, q=c(0.005, 0.995)){
  fits <- list(...)
  if(!all(unlist(lapply(fits, is.srafit))))
    stop("Some arguments passed are not of class srafit")
  if(is.null(names)) names <- paste0('fit',1:length(fits))
  stopifnot(length(names) == length(fits))
  stopifnot(plot %in% c(TRUE, FALSE))
  bscaled <- uscaled <- list()
  pru <- prb <- list()
  for(i in 1:length(fits)){
    ## get posteriors
    uscaled[[i]] <- data.frame(fit=names[i], status=fits[[i]]$uscaled[length(fits[[i]]$year),])
    bscaled[[i]] <- data.frame(fit=names[i], status=fits[[i]]$bscaled[length(fits[[i]]$year),])
    ## Get the prior for B/BMSY
    mean <- fits[[i]]$penalties$deplete.mean
    std <- fits[[i]]$penalties$deplete.cv
    if(!is.null(mean)){
      bb <- seq(qnorm(q[1], mean=mean, sd=std),
                to=qnorm(q[2], mean=mean, sd=std), len=500)
      prb[[i]] <- data.frame(fit=names[i], metric='B/BMSY', status=bb,
                             height=dnorm(bb, mean, std))
    }
    ## Same for U/UMSY
    mean <- fits[[i]]$penalties$harvest.mean
    std <- fits[[i]]$penalties$harvest.sd
    if(!is.null(mean)){
      uu <- seq(qnorm(q[1], mean=mean, sd=std),
                to=qnorm(q[2], mean=mean, sd=std), len=500)
      pru[[i]] <- data.frame(fit=names[i], metric='U/UMSY', status=uu,
                             height=dnorm(uu, mean, std))
    }
  }
  uscaled <- data.frame(metric='U/UMSY', do.call(rbind, uscaled))
  bscaled <- data.frame(metric='B/BMSY', do.call(rbind, bscaled))
  post <- rbind(uscaled, bscaled)
  prior <- rbind(do.call(rbind, prb), do.call(rbind, pru))
  post$fit <- factor(post$fit, levels=names)
  prior$fit <- factor(prior$fit, levels=names)
  prior <- prior[prior$status >= 0,]
  ## Creat ggplot item using the post and prior data sets
  alpha <- 1/length(fits)
  g <- ggplot() +
    geom_histogram(data=post, aes(status, y=..density.., fill=fit),
                   position='identity', alpha=alpha, bins=30) +
    geom_line(data=prior, aes(x=status, y=height, color=fit), lwd=2, alpha=alpha) +
    facet_wrap('metric', scales='free') + geom_vline(xintercept=1, col='red', lty=2)
  if(!is.null(xlim)) g <- g + xlim(xlim)
  ## plot and return
  if(plot) print(g)
  return(invisible(g))
}

#' Plot posterior histograms of MSY reference points and the log posterior
#' density
#'
#' @template plot_args
#' @export
plot_reference <- function(fit){
  stopifnot(is.srafit(fit))
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
  stopifnot(is.srafit(fit))
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
  stopifnot(is.srafit(fit))
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
  stopifnot(is.srafit(fit))
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
  stopifnot(is.srafit(fit))
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
  stopifnot(is.srafit(fit))
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

