### Utility functions for package
#' Plot the prior draws from a fit, showing which crashed the population
#' and which were kept
#' @param fit A returned list from run.SIR
#' @return Nothing. A plot is created.
#'
plot_draws <- function(fit){
  n <- nrow(fit$draws)
  col <- rep('black', len=n)
  col[fit$crashed] <- 'red'
  col[unique(fit$Keepers)] <- 'green'
  ## want the order to be black red green I think
  col <- factor(col, levels=c('black', 'red', 'green'))
  ind <- order(col)
  pairs(fit$draws[ind,], col=col[ind], upper.panel=NULL)
}


#' Read catch data from file
#' @param File to do
#' @export
ReadCatchData <- function(File){
  c="character"; n="numeric"
  D=read.csv(file=File,header=TRUE,colClasses=c(rep(c,6),rep(n,66)),sep=",")
  CatchData=D[,7:72]; colnames(CatchData)=seq(1950,2015)
  DD=list(StockInfo=D[,1:6],CatchData=CatchData)
  return(DD)
}


#' Read in stocks to consider and priors
#' @param File to do
#' @export
ReadPrior <- function(File){
  c="character"; n="numeric"
  D=read.csv(file=File,header=TRUE,colClasses=c(rep(c,8),rep(n,6)),sep=",")
  return(D)
}
