### Utility functions for package

#' Plot time series ribbon for posterior draws
#'
#' @param x Vector of years
#' @param z Matrix of values (years by draws)
#' @param col The base color without transparency.
#' @param alpha.level This is the credible region, e.g. 0.5 will contain
#'   50% of the posterior draws. You can pass multiple levels to get a
#'   shading effect. The transparaceny is proportional to the credible
#'   region (wider regions are lighter).
#' @param alpha.min Minimum level of shading, for fine control of
#'   coloring. Defaults to 0.
#' @param alpha.max See alpha.min. Defaults to 1.
#' @return Nothing, addings filled polygon to current plot.
#'
add.ribbon <- function(x, z, col, alpha.level, alpha.min=0, alpha.max=1){
    ## Pass this a series of years (x) and a matrix of trajectories (z) and it
    ## adds a polygon to the current plot with whose area contains 1-alpha.level
    alpha.level <- sort(alpha.level)
    for(i in 1:length(alpha.level)){
        alpha <- alpha.level[i]
        alpha.col <- alpha.min+alpha*(alpha.max-alpha.min)
        col.poly <- adjustcolor(col, alpha.f = alpha.col)
        quantiles.temp <-
          as.matrix(t(apply(z, 2, quantile, probs=c(alpha/2,1-alpha/2),name=F, na.rm=T)))
        polygon(x=c(x, rev(x)), y=c(quantiles.temp[,1], rev(quantiles.temp[,2])),
                col=col.poly, border=NA)
        lines(x, apply(z, 2, median), lwd=2, col=col)
    }
}


#' Check whether the penalties list is valid.
#'
#' @param penalties A list passed by the user
#'
#' @return Nothing, an error occurs if a problem is found.
check_penalties <- function(penalties){
  if(!is.null(penalties$carry.dist)){
    if(penalties$carry.dist==3){
      if(is.null(penalties$carry.min) | is.null(penalties$carry.max))
        stop("If using carry uniform distribution min and max must be specified")
    }
  }
}


#' Extract prior information for initial biomass, terminal U/UMSY (ustatus)
#' and terminal B/BMSY (bstatus) from a fit, which is useful for plotting
#' and other tasks.
#'
#' @details This function always returns prior values evaluated in natural
#' space, regardless of the distribution used. Thus lognormal priors will
#' be converted to natural space.
#'
#' @param fit An object of class srafit
#' @param metric Which metric to use from 'carry', 'ustatus', 'bstatus'.
#' @param interval Whether to return the prior credible interval or a
#'   sequence of points evaluated at the prior density.
#' @param percentile The interval converage (e.g., .95 is the 95\% credible
#'   interval)
#' @param n.points The number of points if interval is FALSE
#' @return A vector with (lwr, median, upr) if interval is TRUE, and a
#'   data.frame with prior value (x) and prior density (y)
#'
get_prior <- function(fit, metric, interval=TRUE, percentile=.95,
                      n.points=1000){
  pen <- fit$penalties
  metric <- match.arg(metric, choices=c('carry', 'ustatus', 'bstatus'))
  if(metric == 'carry' & pen$carry.dist==3){
    CI <- c(pen$carry.min, (pen$carry.max+pen$carry.min)/2, pen$carry.max)
    xseq <- seq(pen$carry.min, to=pen$carry.max, len=n.points)
    df <- data.frame(x=xseq, y=1/(pen$carry.max-pen$carry.min))
  } else {
    ## Otherwise they are lognormal or normal
    if(metric == 'ustatus'){
      mu <- pen$ustatus.mean; sigma <- pen$ustatus.sd
      dist <- pen$ustatus.dist
    } else if(metric == 'bstatus') {
      mu <- pen$bstatus.mean; sigma <- pen$bstatus.sd
      dist <- pen$bstatus.dist
    } else if(metric == 'carry'){
      mu <- pen$carry.mean; sigma <- pen$carry.sd
      dist <- pen$carry.dist
    } else {
      stop("Invalid metric option")
    }
    if(is.null(dist)) stop("No penalty found for supplied metric")
    ## Symmetic probabilities
    p1 <- (1-percentile)/2
    p2 <- percentile+(1-percentile)/2
    if(dist==2){
      q1 <- qlnorm(p1, mean=mu, sd=sigma)
      q2 <- qlnorm(p2, mean=mu, sd=sigma)
      CI <- c(q1, qlnorm(.5, mean=mu, sd=sigma), q2)
      xseq <- seq(q1, to=q2, len=n.points)
      df <- data.frame(x=xseq, y=dlnorm(xseq, mu, sigma))
    } else {
      q1 <- qnorm(p1, mean=mu, sd=sigma)
      q2 <- qnorm(p2, mean=mu, sd=sigma)
      CI <- c(q1, qnorm(.5, mean=mu, sd=sigma), q2)
      xseq <- seq(q1, to=q2, len=n.points)
      df <- data.frame(x=xseq, y=dnorm(xseq, mu, sigma))
    }
  }
  ## Return
  if(interval)
    return(CI)
  else
    return(df)
}
