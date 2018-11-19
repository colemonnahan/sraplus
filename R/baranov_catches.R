#' Calculate Catches based on Baranov Equation
#'
#' @param log_f log fishing mortality
#' @param m natural mortality
#' @param sel selectivity at age
#' @param b_a biomass at age
#' @param catch observed catches
#' @param use 1 = optimzation, anything else equals catches
#'
#' @return catch or ss
#' @useDynLib sraplus, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
#'
baranov_catches <- function(log_f, m, sel, b_a, catch, use = 1) {

  f <- exp(log_f)

  f_at_a <- f * sel

  catch_hat <-
    sum((f_at_a) / (f_at_a + m) * b_a * (1 - exp(-(f_at_a + m))))

    ss <- (log(catch) - log(catch_hat)) ^ 2

  if (use == 1) {
    out <- ss
  } else{
    out <- catch_hat

  }

  return(out)
}
