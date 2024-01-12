#' Estimate meta Irreproducible Discovery Rate (mIDR).
#'
#' This function extends the Li et al. (2011) copula mixture model, originally
#' implemented  in  idr::est.IDR,  to  accommodate  any number  of  replicates.
#' It  computes  the  local  IDR  for all  pairwise combinations of replicates.
#' Then it computes a "meta" local IDR score using the formula:
#'                     1 - (1 - idr_1)*...*(1 - idr_C(r,2)),
#' where C(r,2) represents the number of all pairwise combinations of scores.
#' Once the meta local IDR is obtained, the expected IDR scores are obtained
#' in the same way as in the traditional IDR procedure.
#'
#' @param x A numeric matrix with rows representing the number of omic features
#' and columns representing the number of sample replicates. The numeric values
#' should be positive and represent significance (not necessarily p-values).
#' @param mu Starting value for the mean of the reproducible component
#' Numeric.
#' @param sigma Starting value for the standard deviation of the reproducible
#' component Numeric.
#' @param rho Starting value for the correlation coefficient of the
#' reproducible component Numeric.
#' @param p Starting value for the proportion of the reproducible component
#' Numeric.
#' @param eps Stopping criterion. Iterations stop when the increment of the
#' log-likelihood is less than "eps" times the log-likelihood.
#' Defaults to 0.001.
#' @param n_threads Number of threads for parallel computing. Numeric. Defaults
#' to 1.
#' @param max.ite Maximum number of iterations. The default is 30.
#' @return Returns a list of two elements:
#' \describe{
#'   \item{\strong{idr}}{A numeric vector of the local meta IDR
#'   for each observation.}
#'   \item{\strong{IDR}}{A numerical vector of the expected meta IDR
#'    for observations that are as irreproducible or more
#'    irreproducible than the given observations.}
#' }
#' @references Q. Li, J. B. Brown, H. Huang, and P. J. Bickel. (2011)
#' Measuring reproducibility of high-throughput experiments.
#' Annals of Applied Statistics, Vol. 5, No. 3, 1752-1779.
#' @import utils
#' @importFrom idr est.IDR
#' @importFrom future plan
#' @import future.apply
#' @export
#' @examples
#' library(eCV)
#' set.seed(42)
#'
#' # Simulate data.
#' out <- simulate_data(scenario = 1, n_reps = 4, n_features = 1e3)
#'
#' # Set initial parameter values.
#' mu <- 2
#' sigma <- 1.3
#' rho <- 0.8
#' p <- 0.7
#'
#' # Get meta local IDR scores.
#' midr_out <- mIDR(x = out$sim_data, mu, sigma, rho, p)
#' \donttest{
#' library(tidyverse)
#' out$sim_data %>%
#'   as.data.frame() %>%
#'   mutate(`Meta idr` = midr_out$idr) %>%
#'   ggplot(aes(x = `Rep 1`, y = `Rep 2`, color = `Meta idr`)) +
#'   geom_point(size = 1) +
#'  scale_color_gradientn(colors=c("#F4364C", "#D5DADD", "#009CA6" ))+ 
#'   theme_classic()
#' }
mIDR <-
  function(x, mu, sigma, rho, p, eps = 0.001, max.ite = 20, n_threads = 1) {
    # Check parameters values are valid.
    if (!is.numeric(x) || any(x <= 0)) {
      stop("Input 'x' must be a numeric matrix with positive values.")
    }
    if (!is.numeric(n_threads) || any(n_threads <= 0)) {
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
      stop("Input 'eps' must be a single numeric value greater than zero.")
    }
    if (!is.numeric(max.ite) || length(max.ite) != 1 || max.ite <= 0) {
      stop("Input 'max.ite' must be a single positive integer value.")
    }

    # Extract number of samples replicates.
    n_reps <- ncol(x)

    # Get all pairwise combinations.
    rep_combs <- utils::combn(n_reps, 2)

    # Get all pairwise IDR scores
    if (n_threads == 1) {
      idr_score_pairs <-
        apply(rep_combs, 2, function(pair) {
          j1 <- pair[1]
          j2 <- pair[2]
          idr::est.IDR(
            x[, c(j1, j2)],
            mu, sigma, rho, p, eps, max.ite
          )$idr
        })
    } else {
      # Set parallel scheduler.
      future::plan("multisession", workers = n_threads)
      message("Using ", n_threads, " workers.")
      idr_score_pairs <-
        future.apply::future_apply(rep_combs, 2, function(pair) {
          j1 <- pair[1]
          j2 <- pair[2]
          idr::est.IDR(
            x[, c(j1, j2)],
            mu, sigma, rho, p, eps, max.ite
          )$idr
        },
        future.seed = TRUE
        )
      # Close multisession buddies.
      future::plan("sequential")
    }

    # Get meta local IDR score.
    idr <- 1 - apply(X = 1 - idr_score_pairs, MARGIN = 1, FUN = prod)

    # Estimate expected meta-IDR.
    o <- order(idr)
    idr.o <- idr[o]
    idr.rank <- rank(idr.o, ties.method = "max")
    IDR.o <- sapply(
      X = idr.rank,
      FUN = function(index, x) {
        mean(x[seq_len(index)])
      },
      x = idr.o
    )
    IDR <- idr
    IDR[o] <- IDR.o
    return(list(
      idr = idr,
      IDR = IDR
    ))
  }
