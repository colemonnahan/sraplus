
#' Quick function to return priors for the Mexican examples.
#' @param N Number of draws
#' @param iStock Stock number (temporary)
#' @return A data.frame that contains N random draws from all priors
#' @export
draw.priors <- function(N, iStock){
  ## initial depletion prior
  ## InitialPrior=rnorm(nrep,Priors$InitialDepletePrior[iStock],Priors$InitialDepleteCV[iStock])
  ## Switched to a lognormal prior on intial depletion so that we're not
  ## generating negative starting values for biomass. Assuming that the
  ## values given are mean and CV and calculating the SD from that.
  InitialPrior <- rlnorm(nrep, log(Priors$InitialDepletePrior[iStock]),
                         sqrt(log(Priors$InitialDepleteCV[iStock]^2+1)))
  ## catch info
  Catch=as.numeric(CatchData[jStock,])
  Cmax=max(Catch,na.rm=TRUE)
  ## carrying capacity prior
  ## correlation between maximum catch and MSY
  Carry=Priors$Kprior[iStock]*Cmax
  Cprior=runif(nrep,min=Carry/2,max=Carry*2)
  ## #######################
  ## priors from FishLife
  ## #######################
  ##------------ taxonomic information for stock -------------##
  ## draw taxonomic information from depletion prior file with stock info
  Taxon <- Priors[iStock,which(colnames(Priors) %in% c("Class", "Order", "Family", "Genus", "Species"))]
  Taxon[which(Taxon=="")] <- "predictive"
  ##------------ find information from FishLife -------------##
  ## find species within FishLife
  sp <- Search_species(Class=Taxon[,"Class"], Order=Taxon[,"Order"], Family=Taxon[,"Family"], Genus=Taxon[,"Genus"], Species=Taxon[,"Species"], ParentChild_gz=Return$ParentChild_gz)$match_taxonomy
  ## find means and covariance matrix
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
