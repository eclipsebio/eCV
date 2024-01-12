#' Assess Reproducibility via Enhanced Coefficient of Variation
#'
#' This function estimates an "enhanced" coefficient of variation (eCV) to 
#' measure the likelihood of an omic feature being reproducible.  The  eCV 
#' is calculated as |SD^2 - Mean^2| / Mean^2, a metric that decreases with
#' noise among replicates and increases with the mean intensity.
#' 
#' Inferences are made based on the probabilities of eCV values originating
#' from the group of reproducible features.  It assumes  that  reproducible
#' features follow a prior Normal distribution with dimension r (number  of
#' replicates). Pseudo replicates are sampled using a Probabilistic Bootstrap,
#' assuming that the global mean vector and variance-covariance matrix across
#' features are close to the prior's hyperparameters. 
#' eCV is computed for each random sample. 
#' The proportion of times the observed eCV is lower than or equal to the eCV
#' from random samples is then taken as the probability of  the  omic feature
#' belonging to the group of reproducible features.
#'
#' @param x A numeric matrix with rows representing omic features and columns 
#' representing sample replicates. Numeric values should be positive and reflect 
#' significance (not necessarily p-values).
#' @param max.ite Number of samples from the null distribution (numeric). 
#' Defaults to 1e4.
#' @param n_threads Number of threads for parallel computing. Numeric. Defaults 
#' to 1.
#' @return Returns a list with two elements:
#' \itemize{
#' \item \emph{\strong{ecv:}} Numeric vector with
#' the estimated eCV values for each omic feature.
#' \item \emph{\strong{post_prob:}} Posterior probability values.
#' }
#' @import stats
#' @importFrom mvtnorm rmvnorm
#' @importFrom future plan
#' @export 
#' @examples
#' library(eCV)
#' set.seed(42)
#' # Simulate data.
#' out <- simulate_data(scenario = 1,n_features=1e3)
#' 
#' # Run eCV
#' ecv_out <- eCV(x = out$sim_data, max.ite = 100)
#' 
#' \donttest{
#' # Plot results.
#' library(tidyverse)
#' 
#' out$sim_data %>% 
#' as.data.frame() %>% 
#' mutate(`eCV Prob` = ecv_out$post_prob) %>%
#' ggplot(aes(x = `Rep 1`, y = `Rep 2`, color = `eCV Prob`)) +
#'  geom_point(size = 1) + 
#'  scale_color_gradientn(colors=c( "#009CA6", "#D5DADD", "#F4364C"))+ 
#'  theme_classic()
#'  }
eCV <- function(x, max.ite = 1e4, n_threads = 1) {
  
  # Check parameters values are valid.
  if (!is.numeric(x) || any(x < 0)) {
    stop("Input 'x' must be a numeric matrix with positive values.")
  }

  if (!is.numeric(n_threads) || any(n_threads <= 0)) {
    stop("Input 'n_threads' must be a numeric matrix with positive values.")
  }
  
  if (!is.numeric(max.ite) || length(max.ite) != 1 || max.ite <= 0) {
    stop("Input 'max.ite' must be a single positive integer value.")
  }
  
  # Get observed eCV.
  omic_feat_means <- rowMeans(x, na.rm = TRUE)
  omic_feat_sd2 <- apply(X = x, MARGIN = 1, FUN = stats::var, na.rm = TRUE)

  ecv <- ifelse(test = omic_feat_means == 0, 
                yes = omic_feat_sd2,  # Avoid Inf when dividing by zero.
                no = abs(omic_feat_sd2 - omic_feat_means^2) / omic_feat_means)
  
  # Get total number of features.
  n_omic_feat <- length(ecv)
  
  # Estimate prior means.
  mu <- colMeans(x,
                 na.rm = TRUE)
  
  # Estimate prior variance-covariance matrix.
  Sigma <- stats::var(x,
               use = "pairwise.complete.obs")
  
  # Get number of probabilistic bootstrap samples indicating if features come
  # from reproducible group.
  if (n_threads == 1) {
  p_num <- 
    rowSums(do.call(what = cbind,
                    args = lapply(X = seq_len(max.ite),
                    FUN =  ecv_null_sampler,
                    n =n_omic_feat,
                    mu = mu,
                    Sigma = Sigma,
                    ecv = ecv) ))
  } else {
    # Set parallel scheduler.
    future::plan("multisession", workers = n_threads)
    message("Using ", n_threads, " workers.")
    p_num <- 
      rowSums(do.call(what = cbind,
                      args = future.apply::future_lapply(X = seq_len(max.ite),
                                                 FUN =  ecv_null_sampler,
                                                 n = n_omic_feat,
                                                 mu = mu,
                                                 Sigma = Sigma,
                                                 ecv = ecv,
                                                 future.seed = TRUE)))
  
    # Close multisession buddies.
    future::plan("sequential")
  }
  
  # Estimate posterior probability of being reproducible.
  post_prob <- 1 - (1 + p_num) / (1 + max.ite)
  
  return(list("eCV" = ecv, "post_prob" = post_prob))
}
