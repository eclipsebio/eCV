#' Estimate IDR with a General Mixture Model
#'
#' This function builds  upon `idr::est.IDR` to extend the Li et al. (2011)
#' copula mixture model to accommodate  an  arbitrary number of replicates.
#' The  term  "General"  in  this context  alludes to  the assumption of a 
#' general multivariate Normal distribution of dimension "n," equal to the
#' number  of  sample  replicates. This  assumption essentially allows the
#' pseudo-likelihood  approach  in  Li et al. (2011) to be extended to any
#' number of replicates.  This  is  achieved  by  modifying the "E" and "M" 
#' steps  of  an  expectation  maximization algorithm to use a multivariate
#' Normal Distribution instead.
#' 
#' @param x Numeric matrix with rows representing the number of omic features 
#' and columns representing the number of sample replicates.
#' The numeric values must be positive and represent significance 
#' (not necessarily p-values).
#' @param mu Starting value for the mean of the reproducible component. Numeric.
#' @param sigma Starting value for the standard deviation of the reproducible 
#' component.
#' @param rho Starting value for the correlation coefficient of the reproducible 
#' component.
#' @param p Starting value for the proportion of the reproducible component.
#' @param eps Stopping criterion. Iterations stop when the increment of the 
#' log-likelihood is less than eps times the log-likelihood. Default is 0.001.
#' @param max.ite Maximum number of iterations. Default is 30.
#' @return Returns a list of three elements:
#' \describe{
#'   \item{\strong{idr}}{A numeric vector of the local IDR 
#'   (Irreproducible Discovery Rate) for each observation (i.e., estimated 
#'   conditional probability for each observation to belong to the 
#'   irreproducible component).}
#'   \item{\strong{IDR}}{A numerical vector of the expected Irreproducible
#'    Discovery Rate for observations that are as irreproducible or more 
#'    irreproducible than the given observations.}
#'   \item{\strong{est_param}}{Estimated parameters: p, rho, mu, sigma.}
#'
#' }
#' @references Q. Li, J. B. Brown, H. Huang, and P. J. Bickel. (2011)
# Measuring reproducibility of high-throughput experiments.
# Annals of Applied Statistics, Vol. 5, No. 3, 1752-1779.
#' @import utils
#' @import stats
#' @importFrom idr get.pseudo.mix
#' @export
#' @examples
#' # 1. Show that gIDR reduces to classical IDR for n=2.
#' 
#' # Load required packages.
#' library(idr)
#' library(eCV)
#' 
#' # Set seed for RNG.
#' set.seed(42)
#' 
#' # Simulate data.
#' out <- simulate_data(scenario = 2, n_features = 1e3)
#' 
#' # Set initial parameter values.
#' mu <- 2
#' sigma <- 1.3
#' rho <- 0.8
#' p <- 0.7
#' 
#' # Compare IDR and gIDR
#' idr.out <- est.IDR(x = out$sim_data, mu, sigma, rho, p)
#' gidr.out <- gIDR(x = out$sim_data, mu, sigma, rho, p)
#' 
#' # Show the results are the same.
#' all.equal(gidr.out$est_param, idr.out$para)
#' 
#' \donttest{
#' 
#' library(tidyverse)
#' # Plot results.
#' out$sim_data %>% as.data.frame() %>% 
#' mutate(idr = gidr.out$idr) %>%
#' ggplot(aes(x=`Rep 1`,y=`Rep 2`,color=idr)) +
#'  geom_point(size=1) + 
#'  scale_color_gradientn(colors=c("#F4364C", "#D5DADD", "#009CA6" ))+ 
#'  theme_classic()
#' 
#' #2. Show gIDR for n=10.
#' 
#' out <- simulate_data(scenario = 1, n_reps = 10, n_features = 1e3)
#' gidr.out <- gIDR(x = out$sim_data, mu, sigma, rho, p)
#' out$sim_data %>% as.data.frame() %>% 
#' mutate(idr = gidr.out$IDR) %>%
#' ggplot(aes(x = `Rep 1`, y = `Rep 2`, color = idr)) +
#'  geom_point(size = 1) + 
#'  scale_color_gradientn(colors=c("#F4364C", "#D5DADD", "#009CA6" ))+ 
#'  theme_classic()
#'}

gIDR <- function (x, mu, sigma, rho, p, eps = 0.001, max.ite = 30)
{
  # Check parameters values are valid.
  if (!is.numeric(x) || any(x <= 0)) {
    stop("Input 'x' must be a numeric matrix with positive values.")
  }
  if (!is.numeric(mu) || length(mu) != 1) {
    stop("Input 'mu' must be a single numeric value.")
  }
  if (!is.numeric(sigma) || length(sigma) != 1) {
    stop("Input 'sigma' must be a single numeric value.")
  }
  if (!is.numeric(rho) || length(rho) != 1) {
    stop("Input 'rho' must be a single numeric value.")
  }
  if (!is.numeric(p) || length(p) != 1) {
    stop("Input 'p' must be a single numeric value.")
  }
  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0) {
    stop("Input 'eps' must be a single numeric value greater than Zero.")
  }
  if (!is.numeric(max.ite) || length(max.ite) != 1 || max.ite <= 0) {
    stop("Input 'max.ite' must be a single positive integer value.")
  }
  
  
  # Get ECDF.
  x.cdf.func <- lapply(X = seq_len(ncol(x)),
                       FUN = function(j) stats::ecdf(x[, j]))
  
  # Re-scale to avoid potential out-of-bound results.
  x.cdf <- lapply(X = seq_len(ncol(x)),
                  FUN = function(j, a) x.cdf.func[[j]](x[, j]) * a,
                  a = nrow(x) / (nrow(x) + 1))
  
  # Store initial parameters values.
  para <- list(mu = mu, sigma = sigma, rho = rho, p = p)
  
  # Initiate outer loop counter.
  j <- 1
  
  # Initiate stopping counter.
  to.run <- TRUE
  
  # Initiate loglihood holders.
  loglik.trace <- loglik.inner.trace <- NULL
  
  # Initiate pseudo data ("Z") values [Eq. 2.5] in Li et al. (2011).
    ## U = n / (n + 1) * CDF (X) (uniform transform of the data "X").
    ## Z = G_inv(U),
    ## where G(Z) = G1(Z) + G0(Z),
    ## G1(Z) = (p / sigma) * MVD((Z - mu) / sigma)
    ## G0(Z) = (1 - p) * MVD(Z)
  Z <- 
    do.call(what = cbind,
            args = lapply(X = seq_len(ncol(x)),
                          FUN = function(j1) 
                            idr::get.pseudo.mix(x = x.cdf[[j1]],
                                                mu = para$mu,
                                                sigma = para$sigma,
                                                rho = para$rho,
                                                p = para$p)))
  
  while (to.run) {
    # Initiate inner loop counter.
    i <- 1
    while (to.run) {
      # "E" step:
       ## Estimate expectation of "K" (a variable indicating if feature comes 
       ## from the Reproducible or the Irreproducible mixture's component) 
       ## given "Z", the pseudo data.
      K <- gIDR.estep(z = Z,
                      mu =  para$mu,
                      sigma = para$sigma,
                      rho = para$rho, 
                      p = para$p)
      
      # "M" step:
        ## Get the maximum likelihood estimates of each parameter given "K".
      para <- gIDR.mstep(Z, K)
      
      # Store loglihood values.
      if (i > 1) l.old <- l.new
      l.new <- loglihood(Z,
                         para$mu,
                         para$sigma,
                         para$rho,
                         para$p)
      loglik.inner.trace[i] <- l.new
      
      # Check inner loop convergence.
      if (i > 1) {
        to.run <- !conv(loglik.inner.trace[i - 1], 
                        loglik.inner.trace[i],
                        eps)
      }
      i <- i + 1
      
      # Check if the EM gets trapped in an oscillation.
      loglik_osc_test(as.vector(loglik.inner.trace))
    }
    
    # Update pseudo-data "Z" with updated parameters values.
    Z <- do.call(what = cbind,
                 args = lapply(X = seq_len(ncol(x)),
                               FUN = function(j1) {
                                 idr::get.pseudo.mix(x = x.cdf[[j1]],
                                                     mu = para$mu,
                                                     sigma=para$sigma,
                                                     rho = para$rho,
                                                     p = para$p)
                              }))
    
    # Store loglihood values.
    if (j > 1) l.old.outer <- l.new.outer
    l.new.outer <- loglihood(Z, 
                             para$mu, 
                             para$sigma, 
                             para$rho, 
                             para$p)
    loglik.trace[j] <- l.new.outer
    
    # Check outer loop convergence.
    if (j == 1)
      to.run <- TRUE
    else {
      if (j > max.ite) to.run <- FALSE
      else to.run <- !conv(l.old.outer, l.new.outer, eps)
    }
    j <- j + 1
    
    # Check if the EM gets trapped by oscillations.
    loglik_osc_test(as.vector(loglik.trace))
  }
  
  # Estimate local IDR.
  idr <- 1 - K
  
  # Estimate expected IDR.
  o <- order(idr)
  idr.o <- idr[o]
  idr.rank <- rank(x = idr.o, ties.method = "max")
  IDR.o <- sapply(X = idr.rank, 
                  FUN = function(index, x) {
                            mean(x[seq_len(index)])
                  },
                  x = idr.o)
  IDR <- idr
  IDR[o] <- IDR.o
  return(list(est_param = list(p = para$p, 
                               rho = para$rho,
                               mu = para$mu,
                               sigma = para$sigma),
              idr = 1 - K, 
              IDR = IDR))
}
