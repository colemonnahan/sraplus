### srafit class definitions

#' Construct an object of class 'srafit'
#' @param fit A list
#' @return A list of class 'srafit'
#' @export
construct_srafit <- function(fit){
  ## should put some checks in here later
  class(fit) <- 'srafit'
  return(fit)
}

#' Test whether an object has class 'srafit'
#' @param x Object to be tested
#' @return TRUE of FALSE whether it is an srafit object.
#' @export
is.srafit <- function(x) inherits(x, "srafit")

#' Print method for class 'srafit'
#' @param fit Object of class srafit
#' @export
summary.srafit <- function(fit){
  ## for now just print it
  print.srafit(fit)
}

#' Print method for class 'srafit'
#' @param fit Object of class srafit
#' @export
print.srafit <- function(fit){
  cat(paste("Fitted object from sraplus for species:\n"))
  cat(fit$Taxon, "\n")
  cat('With', length(unique(fit$Keepers)), 'unique posterior draws\n')
}
