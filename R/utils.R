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

