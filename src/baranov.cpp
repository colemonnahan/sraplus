#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
using namespace Rcpp;
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double baranov(double log_f, double m, NumericVector sel, NumericVector b_a, double catches) {

  double f = exp(log_f);

  NumericVector f_at_a = f * sel;

  double catch_hat =
    sum((f_at_a) / (f_at_a + m) * b_a * (1 - exp(-(f_at_a + m))));


  return pow(log(catch_hat) - log(catches),2);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R

f_guess = 0.2

NatMort = 0.2

num = 1000 * exp(-.2 * (0:19))

Weight <- 1e-6*(num * 0:19)^3

vuln = c(rep(0,10), rep(1,10))


fitted_f <-
  nlminb(
    log(0.2),
    baranov,
    m = NatMort,
    sel = vuln,
    b_a = num * Weight,
    catch = 100000
  )

fitted_f
*/
