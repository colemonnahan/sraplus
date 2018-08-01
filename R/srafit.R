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
  priors <- data.frame(t(sapply(c('carry', 'initial', 'ustatus', 'bstatus'), function(m)
    get_prior(object, metric=m))))
  names(priors) <- c('2.5%', '50%', '97.5%')
  cat('\nPrior intervals:\n')
  print(priors, digits=2)
  cat('\n')
  probs <- c(0.025, .5, .975)
  keep <- object$Keepers
  carry <- quantile(object$draws$Cprior[keep], probs)
  initial <- quantile(object$draws$InitialPrior[keep], probs)
  uscaled <- quantile(object$uscaled[length(object$years),], probs)
  bscaled <- quantile(object$bscaled[length(object$years),], probs)
  posteriors <- rbind(carry, initial, uscaled, bscaled)
  cat('Posterior intervals:\n')
  print(posteriors, digits=2)
  cat('\n')
}

#' Print method for class 'srafit'
#' @param x Object of class srafit
#' @param ... Other arguments to be passed
#' @export
print.srafit <- function(x, ...){
  cat(paste("Fitted object from sraplus for species:\n"))
  cat(x$Taxon, "\n")
  cat('With', length(unique(x$Keepers)), 'unique posterior draws\n')
}
