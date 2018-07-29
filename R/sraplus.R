#' sraplus: Stochastic stock reduction analysis plus biological priors and
#' penalties
#'
#' In development, do not use!
#'
#' @details
#' This package is designed to extend the SRA analysis framework to include
#'   penalties on depletion and priors on life-history parameters derived
#'   from FishLife.
#'
#' @references
#' Thorson, J. T., S. B. Munch, J. M. Cope, and J. Gao. 2017. Predicting
#' life history parameters for all fishes worldwide. Ecological
#' Applications. 27(8):2262-2276. http://onlinelibrary.wiley.com/doi/10.1002/eap.1606/full
#'
#' Walters, C. J., S. J. D. Martell, et
#' al. (2006). "A stochastic approach to stock reduction analysis."
#' Canadian Journal of Fisheries and Aquatic Sciences 63(1): 212-223.
#'
#' @docType package
#' @name sraplus
#' @importFrom stats dnorm rlnorm rnorm runif
#' @importFrom utils read.csv
#' @importFrom grDevices rgb adjustcolor gray
#' @importFrom graphics abline lines pairs plot points axis box legend
#'   mtext par polygon
#' @importFrom stats optimize qnorm median quantile
#' @importFrom utils tail
NULL
