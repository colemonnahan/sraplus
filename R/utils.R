### Utility functions for package


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
