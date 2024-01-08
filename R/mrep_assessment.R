#' Multi-replicate Reproducibility Assessment.
#' 
#' This function wraps the different methods implemented in the  package  eCV to
#' assess reproducibility of omic feature  values coming from two or more sample
#' replicates. The 'method' argument specifies any of the implemented methods: 
#' "IDR", "gIDR", "mIDR", and "eCV".  
#' 
#' The  "IDR"  method  calls  the  traditional  IDR,  as  implemented  in  the
#' package idr (idr::est.IDR).  Regardless  of the number of replicates given to
#' the function, when method="IDR", only the first two are used.
#' Any of the other methods are meant to be used when r >= 2. Both gIDR and mIDR
#' reduce to traditional IDR if r = 2.
#' 
#' @param x A numeric matrix with rows representing the number of omic features
#' and columns representing the number of sample replicates. The numeric values
#' should be positive and represent significance (not necessarily p-values).
#' @param method The name of the method used to assess reproducibility. 
#' Character. Possible values are "IDR", "gIDR", "mIDR", and "eCV". 
#' Defaults to "eCV".
#' @param param List specifying the initial values for the parameters used by 
#' the specified method. If method is any of the IDR variants, param must 
#' be a named list with  "mu", "sigma", "rho", "p", "eps", and "max.ite". If 
#' method = "eCV", param only needs "max.ite".
#' @param n_threads Number of threads for parallel computing. Numeric. Default 
#' to 1. Only used when method is mIDR or eCV.
#' @return A list with two elements
#' \describe{
#'   \item{rep_index}{A numeric vector with values between zero and one, with
#'   smaller values indicating higher reproducibility}
#'   \item{method}{String storing the name of the method used}
#' }
#' @importFrom idr get.pseudo.mix
#' 
#' @export
#' @examples
#' library(eCV)
#' 
#' # Simulate data
#' set.seed(42)
#' out <- simulate_data(scenario = 2, n_reps = 4, n_features = 1e3)
#' 
#' # Define parameters for each method.
#' params <- list(
#'   eCV = list(max.ite = 100),
#'   gIDR = list(
#'      mu = 2,
#'      sigma = 1.3,
#'      rho =0.8,
#'       p = 0.7,
#'     eps = 1e-3,
#'     max.ite = 50),
#'   mIDR = list(
#'      mu = 2,
#'      sigma = 1.3,
#'      rho =0.8,
#'       p = 0.7,
#'     eps = 1e-3,
#'     max.ite = 50))
#' 
#' # Create a list to store results
#' results <- NULL
#' methods <- c("eCV", "gIDR", "mIDR")
#' for (method in methods) {
#'   rep_index <- mrep_assessment(x = out$sim_data,
#'                                method = method,
#'                                param = params[[method]])$rep_index
#'   new_rows <- data.frame(value = rep_index,
#'                          Method = method,
#'                          group = out$sim_params$feature_group)
#'   results <- rbind(results, new_rows)
#' }
#' 
#' # Plot results
#' \donttest{
#' library(tidyverse)
#' results %>% 
#' mutate(group = ifelse(group == 1,"FALSE","TRUE")) %>%
#' ggplot(aes(x=Method, y = value,fill=group)) + 
#' scale_fill_manual(values = c( "#009CA6" , "#F4364C")) + 
#' geom_boxplot() + 
#' theme_classic() + 
#' labs(y="Reproducibility assessment", fill="Reproducible\nfeature")
#' }
mrep_assessment <- function(x, method="eCV", param, n_threads=1) {
  # Check if the 'param' argument is a list
  if (!is.list(param)) {
    stop("The 'param' argument must be a list.")
  }
  # Check if number of threads is right.
  if (!is.numeric(n_threads) || any(n_threads <= 0)) {
    stop("Input 'x' must be a numeric matrix with positive values.")
  }
  # Check if the 'method' is 'eCV' and 'param' contains only 'max.ite'
  if (method == "eCV" && !all(names(param) == "max.ite")) {
    stop("For method 'eCV', 'param' must contain only 'max.ite'.")
  }
  
  # Check if the 'method' is an IDR variant and 'param' contains required elements
  if (method %in% c("IDR", "gIDR", "mIDR") &&
      !all(names(param) %in% c("mu", "sigma", "rho", "p", "eps", "max.ite"))) {
    stop("For IDR variants, 'param' must contain 'mu', 'sigma', 'rho', 'p', 
         'eps', and 'max.ite'.")
  }
  
  # Check if methods allow for parallel compute.
  if (method %in% c("IDR", "gIDR") && n_threads > 1) {
    warning("'n_thread' ignored when 'method' is IDR or gIDR\n")
  }
  
  if (method == "IDR") {
    rep_index <- 
      idr::est.IDR(x = x, 
                   mu = param$mu, 
                   sigma = param$sigma,
                   rho = param$rho,
                   p = param$p,
                   eps = param$eps,
                   max.ite = param$max.ite)$idr
  }
  if (method == "gIDR") {
    rep_index <- 
      gIDR(x = x, 
      mu = param$mu, 
      sigma = param$sigma,
      rho = param$rho,
      p = param$p,
      eps = param$eps,
      max.ite = param$max.ite)$idr
  }
  if (method == "mIDR") {
    rep_index <- 
      mIDR(x = x, 
           mu = param$mu, 
           sigma = param$sigma,
           rho = param$rho,
           p = param$p,
           eps = param$eps,
           max.ite = param$max.ite,
           n_threads = n_threads)$idr
  }
  if (method == "eCV") {
    rep_index <- 
      1 - eCV(x = x, 
           max.ite = param$max.ite,
           n_threads = n_threads)$post_prob
  }
  
  return(list(rep_index = rep_index, 
              method = method))
}

