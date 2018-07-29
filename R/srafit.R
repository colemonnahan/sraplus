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
#' @param object Object of class srafit
#' @param ... Other arguments to be passed
#' @export
summary.srafit <- function(object, ...){
  ## for now just print it
  print.srafit(object, ...)
}

#' Print method for class 'srafit'
#' @param x Object of class srafit
#' @param ... Other arguments to be passed
#' @export
print.srafit <- function(x, ...){
  cat(paste("Objectted object from sraplus for species:\n"))
  cat(x$Taxon, "\n")
  cat('With', length(unique(x$Keepers)), 'unique posterior draws\n')
}
