#' Stochastic SRA to calculate the trajectory given parameters
#' @param Catch To do
#' @param AgeMat Age at maturity
#' @param Steep Steepness
#' @param NatMort Natural mortality
#' @param AgeMax Maximum age, above which there is a plus group.
#' @param Carry Carrying capacity.
#' @param Weight Vector of weights at age, calculated from growth
#'   parameters a and b.
#' @param InitialDeplete Initial level of depletion.
#' @param AgeVulnOffset Optional parameter for specifying age at
#'   vulnerability, defaulting to AgeMat-1. See ?run.SIR for more info.
#' @param Sigma Process error variance, i.e. N(0,Sigma).
#' @param ProcessError A flag for whether to include process error,
#'   defaults to TRUE.
#' @param simulation A list with elements 'NYears' and 'FishMort'. If
#'   supplied the model will run a simulation under that level of U. It it
#'   designed to be used to test equilibrium conditions but can be expanded
#'   later for more expansive simulation testing.
#' @export
AgeModel <- function(Catch,
                     AgeMat,
                     Steep,
                     NatMort,
                     AgeMax,
                     Carry,
                     Weight,
                     InitialDeplete,
                     Sigma,
                     AgeVulnOffset = -1,
                     ProcessError = TRUE,
                     simulation = NULL,
                     carry_form = "ssb") {

  stopifnot(AgeMat > 0)
  AgeVuln <- AgeMat + AgeVulnOffset
  stopifnot(AgeVuln >= 0)

  stopifnot(AgeVuln < AgeMax)
  use.sim <- !is.null(simulation)
  if (use.sim) {
    stopifnot(is.numeric(simulation$FishMort))
    stopifnot(simulation$FishMort >= 0)
    stopifnot(simulation$FishMort <= 1)
    NYears <- simulation$NYears
  } else {
    NYears <- length(Catch)
  }
  ## The process error deviations
  if (ProcessError) {
    devs <- rnorm(NYears, mean = 0, sd = Sigma)
  } else {
    devs <- rep(0, len = NYears)
    Sigma <- 0 # since used below in correction
  }
  ## maximum age
  maxage <- AgeMax
  ## natural survival
  surv  <-  exp(-NatMort)
  ## num is numbers at age; vuln is vulnerability
  num <- vuln <- rep(1, maxage)
  ## 100% selectivity option
  if (AgeMat > 1) {
    ## assume knife edge selectivity
    vuln[1:(AgeVuln)] <- 0
  }
  ## assume logistic maturity with 50% at AgeMat
  mature  <-  1 / (1 + exp(AgeMat - 1:AgeMax))

  rzero_fit <-
    nlminb(
      log(Carry),
      tune_r0,
      carry = Carry,
      m = NatMort,
      mature = mature,
      weight = Weight,
      max_age = AgeMax,
      carry_form = carry_form
    )

  ## per recruit

  rzero <- exp(rzero_fit$par)

  num <- rzero * exp(-NatMort * (0:(maxage - 1)))

  num[maxage] <-  num[maxage - 1] * (exp(-NatMort) / (1 - exp(-NatMort)))

  ## equilibrium spawning biomass
  bzero <- sum(Weight * num)# biomass per recruit

  ssbzero <- sum(Weight * num * mature)# biomass per recruit

  ## include process error
  # rzero  <-  rzero * exp(-Sigma * 2 / 2) #correct for recruitment devs
  ## Beverton-Holt stock-recruit parameters
  alpha <- (ssbzero * (1 - Steep)) / (4 * Steep * rzero)
  betav <- (5 * Steep - 1) / (4 * Steep * rzero)
  ## set initial numbers at age

  r_depletion_fit <- nlminb(
    log(Carry),
    tune_r0,
    carry = Carry * InitialDeplete,
    m = NatMort,
    mature = mature,
    weight = Weight,
    max_age = AgeMax
  )

  r_depletion <- exp(r_depletion_fit$par)

  num <- r_depletion * exp(-NatMort * (0:(maxage - 1)))

  num[maxage] <-  num[maxage - 1] * (exp(-NatMort) / (1 - exp(-NatMort)))

  num0 <- num
  ##  now loop over time
  ## bio is spawning biomass; pop is vulnerable biomass; rec is recruits;
  ## hrstore is the harvest rate
  Vpop <- hrstore <- f_y <- bio <- pop <- rec <- eggs <-
    rep(NA, length = NYears)
  crashed <- FALSE
  for (y in 1:NYears) {
    ## vulnerable biomass
    Vpop[y] <- sum(num * vuln * Weight)
    ## if simulating a population, calculate Catch
    if (use.sim) {
      Catch[y] <- simulation$FishMort * Vpop[y]
    }
    ## harvest rate based on catch

    if (Catch[y] >=  Vpop[y]) {
      break
    }

    f_guess <- -log(1 - Catch[y] / Vpop[y])

    browser()


    baranov_model <-
      TMB::MakeADFun(
        list( m = NatMort,
              sel = vuln,
              b_a = num * Weight,
              catch = Catch[y]),
        list(log_f = log(f_guess)),
        DLL = "baranov")

    fitted_f <- nlminb(baranov_model$par,
                       baranov_model$fn,
                       baranov_model$gr)


    fitted_f <-
      nlminb(
        log(f_guess),
        baranov_catches,
        m = NatMort,
        sel = vuln,
        b_a = num * Weight,
        catch = Catch[y]
      )

    # baranov_catches(fitted_f$par, m = NatMort, sel = vuln, b_a = num * Weight,use = 2, catch = 2)

    if (fitted_f$objective >1) {
      break
    }

    f <- exp(fitted_f$par)

    f_y[y] <- f

    # hr <- Catch[y]/Vpop[y]
    ## hr <- min(.9,Catch[y]/Vpop)

    hr <- 1 - exp(-f)

    hrstore[y] <- hr
    # if(hr>1) break ## break if caught more than available
    ## vulnerable biomass
    pop[y]  <-  sum(num * vuln * Weight)
    if (pop[y] < 0){
      break ## break if negative
    }
    ## spawning biomass
    eggs[y]  <-  sum(num * mature * Weight)
    ## recruits with process error
    rec[y]  <-  eggs[y] / (alpha + betav * eggs[y])  * exp(devs[y])
    last_max <- num[maxage]

    num[2:(maxage)] <- num[1:(maxage - 1)] * exp(-(f * vuln[1:(maxage - 1)] + NatMort))

    num[maxage] <-
      num[maxage] + last_max * exp(-(f * vuln[maxage] + NatMort))

    if (y == 1) {
      num[1] <- num[1]
    } else {
      ## recruitment
      num[1] <-  rec[y - 1]
    }
    if (any(num < 0)){
      browser()
      break

    }
  } #end of loop over time
  ## Calculate UMSY
  ##




  get.equilibrium.catch <- function(U, biomass = FALSE) {
    ## biomass flag is whether to return cmsy or bmsy, see below
    ## equilibrium numbers at age under rate U
    num[1]  <-  rzero
    for (a in 2:(maxage - 1))
      num[a]  <-  (1 - vuln[a - 1] * U) * num[a - 1] * surv
    ## plus group
    num[maxage]  <-
      num[maxage - 1] * ((1 - vuln[maxage - 1] * U) * surv / (1 - (1 - vuln[maxage -
                                                                              1] * U) * surv))
    ## spwaning biomass per recruit fishing at U
    SBPR <- sum(Weight * mature * num) / rzero
    ## recruits in equilibrium fishing at U
    ## R <- max(0, (SBPR-alpha)/(betav*SBPR))
    ## We actually want R to be negative if crashed, so that the function
    ## is continuously differentiable even when crashed. We'll need to
    ## catch this later if CMSY or BMSY is negative.
    R <-  (SBPR - alpha) / (betav * SBPR)
    ## and equlibrium catch per recruit?
    Ca <- sum(num * vuln * Weight * U)
    YPR <- Ca / rzero
    ## Calculate catch for all recruits
    if (biomass)
      ## equilibrium biomass is vulnerable biomass per recruit times
      ## equilibrium recruits
      return (sum(Weight * num * vuln) / rzero * R)
    else
      ## equilibrium catch is yield per recruit times equilibrium recruits
      return(R * YPR)
  }
  ## If a realistic trajectory calculate MSY
  if (any(is.na(pop))) {
    crashed <- TRUE
    bmsy <- NA
    fit <- list(par = NA, objective = NA)
  } else {

    # fit2 <-
    #   optimize(
    #     get.equilibrium.catch,
    #     interval = c(0, 1),
    #     maximum = TRUE,
    #     tol = .001
    #   )

    fit <-
      nlminb(
        log(NatMort),
        get_eq_catch,
        rzero = rzero,
        NatMort = NatMort,
        maxage = maxage,
        vuln = vuln,
        Weight = Weight,
        mature = mature,
        alpha = alpha,
        betav = betav
      )

    ## after finding cmsy put it back in to get bmsy
    # bmsy <- get.equilibrium.catch(fit$maximum, biomass = TRUE)

    bmsy <-get_eq_catch(
        fit$par,
        rzero = rzero,
        NatMort = NatMort,
        maxage = maxage,
        vuln = vuln,
        Weight = Weight,
        mature = mature,
        alpha = alpha,
        betav = betav,
        biomass = TRUE
      )

  }
  out <-
    list(
      pop = pop,
      Vpop = Vpop,
      hr = hrstore,
      f_y = f_y,
      umsy = 1 - exp(-exp(fit$par)),
      fmsy = exp(fit$par),
      cmsy = -fit$objective,
      bmsy = bmsy,
      crashed = crashed,
      devs = devs
    )
  if (use.sim) {
    u.seq <- seq(0, 1, len = 100)
    c.seq <-
      sapply(u.seq, function(u)
        max(0, get.equilibrium.catch(u)))
    b.seq <-
      sapply(u.seq, function(u)
        max(0, get.equilibrium.catch(u, biomass = TRUE)))
    o <- list(
      num0 = num0,
      num = num,
      Catch = Catch,
      Vpop = Vpop,
      CatchEquilibrium = get.equilibrium.catch(simulation$FishMort),
      BiomassEquilibrium = get.equilibrium.catch(simulation$FishMort, TRUE),
      u.seq = u.seq,
      c.seq = c.seq,
      b.seq = b.seq
    )
    out <- c(out, o)
  }
  return(out)
} #end of function
