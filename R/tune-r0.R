tune_r0 <-
  function(log_r0,
           carry,
           m,
           max_age,
           weight,
           mature,
           sigma = 0.2,
           kform = "ssb",
           rcorrect = FALSE) {

    if (rcorrect == T){

      r0 <- exp(log_ro - sigma^2/2)
    }

    r0 <- exp(log_r0)

    n0 <- r0 * exp(-(0:(max_age - 1)) * m)

    n0[max_age] <- n0[max_age - 1] * (exp(-m) / (1 - exp(-m)))

    if (kform == "ssb") {
      b0 <- sum(n0 * weight * mature)
    } else{
      b0 <- sum(n0 * weight)

    }

    out <- (carry - b0) ^ 2

  }
