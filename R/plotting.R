### Plotting functions for package

#' Plot priors vs posteriors for carrying capacity, initial depletion,
#' terminal year U/UMSY and B/BMSY for a series of fits.
#' @template plot_args
#' @param percentile Optional percentile to calculate prior
#'   density range. Defaults to 99\%.
#' @return An invisible ggplot2 object which can be printed or saved.
#' @export
plot_penalties <- function(..., names=NULL, plot=TRUE,
                           percentile=0.99){
  fits <- list(...)
  if(!all(unlist(lapply(fits, is.srafit))))
    stop("Some arguments passed are not of class srafit")
  if(is.null(names)) names <- paste0('fit',1:length(fits))
  stopifnot(length(names) == length(fits))
  stopifnot(plot %in% c(TRUE, FALSE))
  poc <- poi <- pou <- pob <- list()
  prc <- pri <- pru <- prb <- list()
  for(i in 1:length(fits)){
    ## get posteriors
    keep <- fits[[i]]$Keepers
    poc[[i]] <- data.frame(fit=names[i], metric='carry', posterior=fits[[i]]$draws$Cprior[keep])
    poi[[i]] <- data.frame(fit=names[i], metric='initial', posterior=fits[[i]]$draws$InitialPrior[keep])
    pou[[i]] <- data.frame(fit=names[i], metric='ustatus', posterior=fits[[i]]$uscaled[length(fits[[i]]$year),])
    pob[[i]] <- data.frame(fit=names[i], metric='bstatus', posterior=fits[[i]]$bscaled[length(fits[[i]]$year),])
    ## Get the prior for carry
    if(!is.null(fits[[i]]$penalties$carry.dist)){
      cc <- get_prior(fits[[i]], metric='carry', interval=FALSE, percentile=percentile)
      prc[[i]] <- data.frame(fit=names[i], metric='carry', prior=cc[,1], density=cc[,2])
    }
    ## initial
    if(!is.null(fits[[i]]$penalties$initial.dist)){
      ii <- get_prior(fits[[i]], metric='initial', interval=FALSE, percentile=percentile)
      pri[[i]] <- data.frame(fit=names[i], metric='initial', prior=ii[,1], density=ii[,2])
    }
    ## Same for U/UMSY
    if(!is.null(fits[[i]]$penalties$ustatus.dist)){
      uu <- get_prior(fits[[i]], metric='ustatus', interval=FALSE, percentile=percentile)
      pru[[i]] <- data.frame(fit=names[i], metric='ustatus', prior=uu[,1], density=uu[,2])
    }
    ## B/BMSY
    if(!is.null(fits[[i]]$penalties$bstatus.dist)){
      bb <- get_prior(fits[[i]], metric='bstatus', interval=FALSE, percentile=percentile)
      prb[[i]] <- data.frame(fit=names[i], metric='bstatus', prior=bb[,1], density=bb[,2])
    }
  } ## end loop over fits
  ## Massage these into ggplot format
  post <- do.call(rbind, c(poc, poi, pou, pob))
  prior <- do.call(rbind, c(prc, pri, prb, pru))
  post$fit <- factor(post$fit, levels=names)
  prior$fit <- factor(prior$fit, levels=names)
  ## None of these should be negative so chop it off
  prior <- prior[prior$prior >= 0,]
  ## Creat ggplot item using the post and prior data sets
  alpha <- 1/length(fits)
  g <- ggplot() +
    geom_histogram(data=post, aes(posterior, y=..density.., fill=fit),
                   position='identity', alpha=alpha, bins=30) +
    geom_line(data=prior, aes(x=prior, y=density, color=fit), lwd=2, alpha=alpha) +
    facet_wrap('metric', scales='free') + geom_vline(xintercept=1, col='red', lty=2)
  ## plot and return
  if(plot) print(g)
  return(invisible(g))
}

#' Plot posterior histograms of MSY reference points and the log posterior
#' density
#'
#' @template plot_args
#' @export
plot_reference <- function(..., names=NULL, plot=TRUE){
  fits <- list(...)
  if(!all(unlist(lapply(fits, is.srafit))))
    stop("Some arguments passed are not of class srafit")
  if(is.null(names)) names <- paste0('fit',1:length(fits))
  stopifnot(length(names) == length(fits))
  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))
  tmp <- list()
  for(i in 1:length(fits)){
    x <- fits[[i]]
    tmp[[i]] <- data.frame(fit=names[i], BMSY=x$bmsy, UMSY=x$umsy, MSY=x$cmsy,
                      log_likelihood=log(x$likes[x$Keepers]))
  }
  out <- do.call(rbind, tmp)
  results.long <- reshape2::melt(out, 'fit', variable.name='metric')
  g <- ggplot(results.long) +
    geom_histogram(aes(value, y=..density.., fill=fit),
                   position='identity', alpha=1/length(fits), bins=30) +
    facet_wrap('metric', scales='free')
  if(plot) print(g)
  return(invisible(g))
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
  ci <- get_prior(fit, metric='bstatus', interval=TRUE)
  lines(x=rep(tail(fit$year,1),2), ci[c(1,3)], col=2, lwd=2)
  points(x=tail(fit$year,1), y=ci[2], col=2, cex=1.5, pch=16)
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
  ci <- get_prior(fit, metric='ustatus', interval=TRUE)
  lines(x=rep(tail(fit$year,1),2), ci[c(1,3)], col=2, lwd=2)
  points(x=tail(fit$year,1), y=ci[2], col=2, cex=1.5, pch=16)
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
  ## Add confidence interval
  ci <- get_prior(fit, metric='carry', interval=TRUE)
  lines(x=rep(head(fit$year,1),2), ci[c(1,3)], col=2, lwd=2)
  points(x=head(fit$year,1), y=ci[2], col=2, cex=1.5, pch=16)
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


#' Plot management quantities for one or more fits
#' @template plot_args
#' @param lims The limits for metrics. A single range can be passed and
#'   used for all 4 metrics, or a list of ranges. NULL specifies to
#'   determine the ranges for each metric from the data.
#' @export
plot_fit <- function(..., names=NULL, lims=c(0,3)){
  axis.col <- gray(.5)
  cex.label <- .7
  box.tmp <- function() box(col=axis.col)
  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))
  fits <- list(...)
  if(!all(unlist(lapply(fits, is.srafit))))
    stop("Some arguments passed are not of class srafit")
  if(is.null(names))
    names <- paste0('fit',1:length(fits))
  n <- 4
  nfits <- length(fits)
  cols <- 1:nfits
  par(mfcol=c(n,n), mar=c(.1,.1,.1,.1), yaxs="i", xaxs="i", mgp=c(.25, .25,0),
      tck=-.02, cex.axis=.65, col.axis=axis.col, oma=c(2.5, 2.5, .5,.5))
  label.cex <- .8
  ## Get quantities that span all fits
  q1 <- .95
  if(is.null(lims)){
  ## ylims.ts <- list(
  ##   c(0,quantile(unlist(lapply(fits, function(x) x$depletion)), probs=q1)),
  ##   c(0,quantile(unlist(lapply(fits, function(x) x$bscaled)), probs=q1)),
  ##   c(0,quantile(unlist(lapply(fits, function(x) x$U)), probs=q1)),
  ##   c(0,quantile(unlist(lapply(fits, function(x) x$uscaled)), probs=q1)))
  lims <- list(
    1.1*c(0,max(unlist(lapply(fits, function(x) x$depletion[nrow(x$depletion),])))),
    1.1*c(0,max(unlist(lapply(fits, function(x) x$bscaled[nrow(x$depletion),])))),
    1.1*c(0,max(unlist(lapply(fits, function(x) x$U[nrow(x$depletion),])))),
    1.1*c(0,max(unlist(lapply(fits, function(x) x$uscaled[nrow(x$depletion),])))))
  } else {
    ## If single range provided, replicate it for each metric
    if(!is.list(lims)){
      lims <- rep(list(lims), n)
    } else {
      ## Otherwise user passed list so check it
      stopifnot(length(lims) == n)
    }
  }
  for(cc in 1:n){
    for(rr in 1:n){
      for(ii in 1:length(fits)){
        x <- fits[[ii]]
        term <- cbind(depletion=x$depletion[nrow(x$depletion),], bscaled=x$bscaled[nrow(x$bscaled),],
                      U=x$U[nrow(x$U),],
                      uscaled=x$uscaled[nrow(x$uscaled),])
        time <- list(depletion=x$depletion, bscaled=x$bscaled, U=x$U, uscaled=x$uscaled)
        term.names <- c("Depletion", "B/BMSY", "U", "U/UMSY")
        ## case 1 is diagonal; time series
        tmp <- time[[rr]]
        if(rr==cc){
          if(ii==1){
            plot(x$years, tmp[,1], type='n', col=cols[ii],
                 ylim=lims[[rr]], axes=FALSE, ann=FALSE)
            }
            mtext(text=term.names[rr], line=-1.5, cex=cex.label)
            add.ribbon(x=x$years, z=t(time[[cc]]), col=cols[ii], alpha.level=c(.5))
            abline(h=1, lty=3, col=axis.col)
            box.tmp()
        }
        ## case 2 is lower off diagonal: scatter plot of terminal years
        if(rr > cc) {
          if(ii==1){
            plot(term[,cc], term[,rr],  col=cols[ii], ylim=lims[[rr]],
                 xlim=lims[[cc]], ann=FALSE,
                 xlab=term.names[cc], ylab=term.names[rr], axes=FALSE)
            abline(v=1, h=1, col=axis.col, lty=3)
            box.tmp()
          } else {
            points(term[,cc], term[,rr],  col=cols[ii])
          }
        }
        ## case 3 is upper diagonal and is blank for now except for some
        ## special cases
        if(rr < cc) {
          if(ii==1){
            plot(x$years, x$years, ylim=c(0, 1.1*max(x$Catch)), type='n', ann=FALSE, axes=FALSE)
          }
          ## Catch in the upper right
          if(cc==n & rr ==1){
            lines(x$years, x$Catch)
            axis(1); axis(2)
            box.tmp()
            mtext('Catches', side=2, line=1.5, cex=cex.label)
          }
          if(cc==n & rr==2){
            tmp <- paste(names,unlist(lapply(fits, function(x)
              length(unique(x$Keepers)))), sep='=')
            legend('center', title='Number of draws', pch=1, col=cols,
                   legend=tmp, bty='n')
          }
        }
      } # end loop over fits
      ## Add special cases of axes on the ends
      if(cc==1)
        mtext(term.names[rr], line=1+ifelse(rr %% 2 ==1, .1, .1),
              cex=cex.label, side=2)
      if(rr==n) {
        par( mgp=c(.05, ifelse(cc %% 2 ==0, 0, .5),0) )
        axis(1, col=axis.col, lwd=.5)
      }
      if(cc==1 & rr >1) {
        ## par( mgp=c(.05, ifelse(rr %% 2 ==1, .15, .65),0) )
        axis(2, col=axis.col, lwd=.5)
      }
      if(rr==1 & cc ==1){
        ##    par( mgp=c(.05, ifelse(rr %% 2 ==1, .15, .65),0) )
        axis(2, col=axis.col, lwd=.5)
      }
      ## if(rr==n)
      ##   mtext(term.names[cc], side=1, line=1.5+ifelse(rr %% 2 ==1, 0, .01),
      ##         cex=cex.label)
    } # end of loop over rows
  } # end of loop over columns
}

