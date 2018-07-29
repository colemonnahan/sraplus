
#' Draw samples from the prior distribution as defined by biological
#'   parameters from FishLife and initial depletion distribution.
#' @param N Number of draws
#' @param InitialDepletePrior Mean of initial lognormal depletion penalty
#' @param InitialDepleteCV CV of intial lognormal depletion penalty
#' @param Kprior Maybe a multiplier on K?
#' @param Catch Vector of catch
#' @param Taxon A list with taxonomic information, including 'Class',
#'   'Order', 'Family', 'Genus', 'Species'. Any blanks will be converted to
#'   'predictive' internally when querying FishLife.
#' @param Kscale Scalar to control the initial biomass which is generated
#'   as runif(N, Carry/Kscale, Carry*Kscale) where Carry=Kprior*max(Catch).
#' @return A list that contains 'draws', a data.frame with N random draws from all
#'   priors, and 'Taxon' which is a character string that is returned from
#'   the FishLife query match.
draw.priors <- function(N, InitialDepletePrior, InitialDepleteCV, Kprior,
                        Catch, Taxon, Kscale=2){
  ## initial depletion prior
  ## Switched to a lognormal prior on intial depletion so that we're not
  ## generating negative starting values for biomass. Assuming that the
  ## values given are mean and CV and calculating the SD from that.
  InitialPrior <-
    rlnorm(N, log(InitialDepletePrior), sqrt(log(InitialDepleteCV^2+1)))
  ## catch info
  Cmax <- max(Catch,na.rm=TRUE)
  ## carrying capacity prior
  ## correlation between maximum catch and MSY
  Carry <- Kprior*Cmax
  Cprior <- runif(N,min=Carry/Kscale,max=Carry*Kscale)
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
  Which <- grep(sp[1], Return$ParentChild_gz[,"ChildName"])
  Mean <- Return$beta_gv[Which,]
  Cov <- Return$Cov_gvv[Which,,]
  ##------------ deviates from multivariate normal -------------##
  ## parameters for multvariate normal distribution
  ##  params_mvn <- c("tm", "M", "tmax", "Loo", "K", "Winfinity") #"ln_var", "Loo", "h", "logitbound_h")
  params_mvn <- c("tm", "M", "tmax", "Loo", "K", "Winfinity","ln_var", "logitbound_h")
  Mean_mvn <- Mean[which(names(Mean) %in% params_mvn)]
  Cov_mvn <- Cov[which(rownames(Cov) %in% params_mvn), which(colnames(Cov) %in% params_mvn)]
  ## draw random deviates from multivariate normal distribution between
  draws_mvn <- data.frame(rmvnorm(N, mean=Mean_mvn, sigma=Cov_mvn))
  ## exponentiate and add length-weight params
  draws_mvn <- draws_mvn %>%
    mutate(Loo = exp(Loo)) %>%
    mutate(K = exp(K)) %>%
    mutate(Winfinity = exp(Winfinity)) %>%
    mutate(tmax = exp(tmax)) %>%
    mutate(tm = exp(tm)) %>%
    mutate(M = exp(M)) %>%
    mutate("lwb" = 3.04)  %>% ## from Froese, Thorson, and Reyes meta-analysis, mean value for all fish
    mutate("lwa" = Winfinity / (Loo ^ lwb)) %>%
    mutate("h"=  0.2+ 0.8*(1/(1+exp(-logitbound_h)))) %>%
    mutate("Sigma"=exp(ln_var)) %>%
    mutate("Cprior"=Cprior) %>%
    mutate("InitialPrior"=InitialPrior)
  draws <- draws_mvn
  draws$logitbound_h <- draws$ln_var <- NULL
  if(any(draws$Sigma > 1.5)){
    warning(paste(sum(draws$Sigma > 1.5),
                  "of Sigma values were truncated to 1.5"))
    draws$Sigma[draws$Sigma > 1.5] <-  1.5
  }
  if(any(draws$M > 5)){
    warning(paste(sum(draws$M > 5),
                  "of M values were truncated to 3"))
    draws$M[draws$M > 5] <-  5
  }
  return(list(draws=draws, Taxon=sp[1]))
}
