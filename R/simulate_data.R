#' Simulates omic features into reproducible and irreproducible groups
#' 
#' This function  is an extension  of the copula mixture  model  simulations 
#' presented in Li et al. (2011). It  generates samples of n_features  pairs
#' of omic features for n_reps (>=2)  replicates.  The  state  of each  omic
#' feature  (i.e., reproducible or irreproducible) is determined by sampling
#' from a binomial variable K with a vector of probabilities, P.  
#' The  vector  P  represents the mixing probability between two multivariate
#' normal distributions. The elements of P are associated with reproducibility. 
#' For example, if  K  can only assume two values, say 0 or 1,  then  K  can 
#' represent groups of reproducible or irreproducible features.
#' 
#' The dimension  of each normal distribution is determined by the number of
#' replicates,  r.  The  "scenario"  argument  controls  the  values  of the
#' parameters according  to  the  simulation  scenarios outlined in Li et al.
#' (2011) (Table I in the article).  Scenario  1  corresponds to a situation 
#' where reproducible features are highly correlated  and  exceed the number
#' of irreproducible features.  Scenario 2  corresponds to a situation where
#' the  reproducible  features  are  less  than  the irreproducible ones and
#' exhibit   low   correlation.   Scenario 3   represents  situations  where
#' reproducible  features are less than irreproducible ones but still highly
#' correlated. Scenario 4  is  a  generalization  of  Scenario 1,  with  the
#' addition  of  a  component  of  “reproducible noise” consisting of highly
#'  correlated but low-intensity features.
#' 
#' @param n_reps Number of sample replicates. Numeric. Defaults to 2.
#' @param n_features Number of omic features to simulate. Numeric. 
#' Defaults to 1e4.
#' @param scenario Combination of parameters' values defining scenarios in 
#' Li et al. (2011). Numeric. Possible values are 1, 2, 3, or 4. 
#' Defaults to 1.
#' @return Returns a list of two elements:
#' \itemize{
#' \item \emph{\strong{sim_data:}} Matrix of dimensions n_features x n_reps with
#' the simulated numerical values for each feature.
#' \item \emph{\strong{sim_params:}} List with all the parameter values.
#' }
#' @references Q. Li, J. B. Brown, H. Huang, and P. J. Bickel. (2011)
#'  Measuring reproducibility of high-throughput experiments. 
#'  Annals of Applied Statistics, Vol. 5, No. 3, 1752-1779.
#' @import stats
#' @importFrom mvtnorm rmvnorm
#' @export
#' @examples
#' library(eCV)
#' set.seed(42)
#' out <- simulate_data(scenario = 1)
#' \donttest{
#' library(tidyverse)
#' out$sim_data %>% as.data.frame() %>% 
#' mutate(`Features group` = as.character(out$sim_params$feature_group)) %>%
#' ggplot(aes(x=`Rep 1`,y=`Rep 2`,color=`Features group`)) +
#'  geom_point(size=1, alpha=0.5) + 
#'  scale_color_manual(values = c( "#009CA6" , "#F4364C")) + 
#'  theme_classic()
#' }
simulate_data <- function(n_reps = 2,
                          n_features = 10000, 
                          scenario = 1) {
  
  # Check parameter values.
  if (!is.numeric(scenario) && !scenario %in% 1:4) 
    stop("scenario must be an integer value from 1 to 4.")

  if (!is.numeric(n_reps) || n_reps < 2) 
    stop("n_reps must be an integer greater than or equal to 2.")

  if (!is.numeric(n_features) || n_features < 1) 
    stop("n_features must be a positive integer.")

  # Set parameter values according to the scenario used.
  if (scenario == 1) {
    p <- c(0.35, 0.65)
    mu <- c(0, 2.5)
    rho <- c(0, 0.84)
  } 
  if (scenario == 2) {
    p <- c(0.7, 0.3)
    mu <- c(0, 2.5)
    rho <- c(0, 0.4)
  }
  if (scenario == 3) {
    p <- c(0.95, 0.05)
    mu <- c(0, 2.5)
    rho <- c(0, 0.84)
  }
  if (scenario == 4) {
    p <- c(0.28, 0.65, 0.07)
    mu <- c(0, 3, 0)
    rho <- c(0, 0.4, 0.64)
  }
  sigma2 <- rep(1, length(p))
  
  # Simulate features.
  n_features_by_group <- round(n_features * p)
  features_groups <- rep(seq_along(p), n_features_by_group)
  sim_features <- 
    do.call(what = rbind, 
            args = lapply(X = seq_along(p), 
                          FUN = function(k) {
                              # Create vector of means.
                              Mu <- rep(mu[k], n_reps)
                              # Create variance-covariance matrix.
                              Sigma <- matrix(data = rho[k], 
                                              nrow = n_reps, 
                                              ncol = n_reps)
                              diag(Sigma) <- rep(x = sigma2[k],
                                                 n_reps)
                              # Sample features values from a MVN distribution.
                              Z <- mvtnorm::rmvnorm(n = n_features_by_group[k],
                                                    mean = Mu,
                                                    sigma = Sigma)
                              # Transform to positive values.
                              -log(1/ (1 + exp(Z)))
    }))
  colnames(sim_features) <- paste("Rep",seq_len(n_reps))
  sim_features <- sim_features[seq_len(n_features), ]
  features_groups <- features_groups[seq_len(n_features)]
  
  # Create a named vector to store parameter values
  sim_params <- list("n_reps" = n_reps, 
                     "n_features" = n_features, 
                     "mu" = mu,
                     "rho" = rho,
                     "sigma2" = sigma2,
                     "feature_group" = features_groups,
                     "p" = p)
  
  # Return simulations and parameters settings.
  return(list(sim_data = sim_features, 
              sim_params = sim_params))
}
