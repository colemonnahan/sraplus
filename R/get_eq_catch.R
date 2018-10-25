#' get equilibrium catch
#'
#' @param log_f log fishing mortality
#' @param rzeroC corrected recruitment deviate
#' @param NatMort natural mortality
#' @param maxage max age
#' @param vuln fishing vulnerability at age
#' @param Weight weight at age
#' @param alpha recruitment alpha
#' @param betav recruitment beta
#' @param biomass logical TRUE = Bmsy, false = catch

#'
#' @return either biomass or catch
#' @export
#'
get_eq_catch <-
  function(log_f,
           rzeroC,
           NatMort,
           maxage,
           vuln,
           Weight,
           mature,
           alpha,
           betav,
           biomass = FALSE) {
    ## biomass flag is whether to return cmsy or bmsy, see below
    ## equilibrium numbers at age under rate U

    f <- exp(log_f)

    num <- rep(NA, maxage)

    num[1]  <-  rzeroC
    for (a in 2:(maxage - 1)) {
      num[a]  <-  num[a - 1] * exp(-(f * vuln[a - 1] + NatMort))
    }

    num[maxage] <-
      (num[maxage - 1] * exp(-(f * vuln[maxage] + NatMort))) / (1 -  exp(-(f * vuln[maxage] + NatMort)))

    ## spwaning biomass per recruit fishing at U
    SBPR <- sum(Weight * mature * num) / rzeroC
    ## recruits in equilibrium fishing at U
    ## R <- max(0, (SBPR-alpha)/(betav*SBPR))
    ## We actually want R to be negative if crashed, so that the function
    ## is continuously differentiable even when crashed. We'll need to
    ## catch this later if CMSY or BMSY is negative.
    R <-  (SBPR - alpha) / (betav * SBPR)
    ## and equlibrium catch per recruit?

    Ca <-
      sum((f * vuln) / (f * vuln + NatMort) * (num * Weight) * (1 - exp(-(f * vuln + NatMort))))

    # Ca <- sum(num * vuln * Weight * U)
    YPR <- Ca / rzeroC
    ## Calculate catch for all recruits
    if (biomass)
      ## equilibrium biomass is vulnerable biomass per recruit times
      ## equilibrium recruits
      return (sum(Weight * num * vuln) / rzeroC * R)
    else
      ## equilibrium catch is yield per recruit times equilibrium recruits
      return(-(R * YPR))
  }
