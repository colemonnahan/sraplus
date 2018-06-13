
#' Quick function to return priors for the Mexican examples.
#' @param N Number of draws
#' @param InitialDepletePrior Mean of initial lognormal depletion penalty
#' @param InitialDepleteCV CV of intial lognormal depletion penalty
#' @param Kprior Maybe a multiplier on K?
#' @param iStock Stock number (temporary)
#' @param Catch Vector of catch
#' @param Taxon A named vector with elements "Class, Order, Family, Genus,
#'   Species". Either provide a name or "" to specify predictive.
#' @return A data.frame that contains N random draws from all priors
#' @export
draw.priors <- function(N, InitialDepletePrior, InitialDepleteCV, Kprior,
                        Catch, Taxon){
  ## initial depletion prior
  ## Switched to a lognormal prior on intial depletion so that we're not
  ## generating negative starting values for biomass. Assuming that the
  ## values given are mean and CV and calculating the SD from that.
  InitialPrior <-
    rlnorm(nrep, log(InitialDepletePrior), sqrt(log(InitialDepleteCV^2+1)))
  ## catch info
  Cmax <- max(Catch,na.rm=TRUE)
  ## carrying capacity prior
  ## correlation between maximum catch and MSY
  Carry <- Kprior*Cmax
  Cprior <- runif(nrep,min=Carry/2,max=Carry*2)
  ## #######################
  ## priors from FishLife
  ## #######################
  ##------------ taxonomic information for stock -------------##
  ## draw taxonomic information from depletion prior file with stock info
  Taxon[which(Taxon=="")] <- "predictive"
  ##------------ find information from FishLife -------------##
  ## find species within FishLife
  sp <- Search_species(Class=Taxon["Class"], Order=Taxon["Order"],
                       Family=Taxon["Family"], Genus=Taxon["Genus"],
                       Species=Taxon["Species"],
                       ParentChild_gz=Return$ParentChild_gz)$match_taxonomy
  ## find means and covariance matrix
  ### Is Return put on global space by Search_species??
  Which <- grep(sp[1], Return$ParentChild_gz[,"ChildName"])
  Mean <- Return$beta_gv[Which,]
  Cov <- Return$Cov_gvv[Which,,]
  ##------------ deviates from multivariate normal -------------##
  ## parameters for multvariate normal distribution
  params_mvn <- c("tm", "M", "tmax", "Loo", "K", "Winfinity") #"ln_var", "Loo", "h", "logitbound_h")
  Mean_mvn <- Mean[which(names(Mean) %in% params_mvn)]
  Cov_mvn <- Cov[which(rownames(Cov) %in% params_mvn), which(colnames(Cov) %in% params_mvn)]
  ## draw random deviates from multivariate normal distribution between
  draws_mvn <- data.frame(rmvnorm(nrep, mean=Mean_mvn, sigma=Cov_mvn))
  ## exponentiate and add length-weight params
  draws_mvn <- draws_mvn %>%
    mutate(Loo = exp(Loo)) %>%
    mutate(K = exp(K)) %>%
    mutate(Winfinity = exp(Winfinity)) %>%
    mutate(tmax = exp(tmax)) %>%
    mutate(tm = exp(tm)) %>%
    mutate(M = exp(M)) %>%
    mutate("lwb" = 3.04)  %>% ## from Froese, Thorson, and Reyes meta-analysis, mean value for all fish
    mutate("lwa" = Winfinity / (Loo ^ lwb))
  ##------------ recruitment deviations -------------##
  Sigma <- rlnorm(nrep, Mean["ln_var"], 0.3)
  Steep <- rlnorm(nrep, log(Mean["h"]), 0.1)
  draws <- draws_mvn %>%
    mutate("Sigma" = Sigma) %>%
    mutate("h" = Steep) %>%
    mutate("Cprior"=Cprior) %>%
    mutate("InitialPrior"=InitialPrior)
  return(draws)
}
