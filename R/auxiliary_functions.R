
# Density of pseudo data ("Z") under Multivariate Normal Distribution (MVD).
## U = n / (n + 1) * CDF (X) (uniform transform of the data "X").
## Z = G_inv(U),
## where G(Z) = G1(Z) + G0(Z),
## G1(Z) = (p / sigma) * MVD((Z - mu) / sigma)
## G0(Z) = (1 - p) * MVD(Z)
norm_dens <- function(z, mu, sigma, rho) 
  {
  # Use SD and correlation to define variance-covariance matrix.
  Sigma <- matrix(rho, ncol(z), ncol(z)) * sigma^2
  diag(Sigma) <- sigma^2 

  mvtnorm::dmvnorm(z, 
                   log = TRUE,
                   checkSymmetry = FALSE,
                   mean = rep(mu, ncol(z)),
                   sigma = Sigma)
}

# Predict assignment to reproducible component ("K") given parameters values.
## Expectation step.
gIDR.estep <- function (z, mu, sigma, rho, p)
{
  p / ((1 - p) * exp(norm_dens(z, mu = 0, sigma = 1, rho = 0) - 
                            norm_dens(z, mu, sigma, rho)) + p)
}

# Estimate parameters values via Maximum Likelihood given "K".
## Maximization step.
gIDR.mstep <- function (z, K)
{
  # Get mixing probability 
  p <- mean(K)
  
  # Get mean of reproducible component.
  mu <- sum(rowMeans(z * K)) / sum(K)
  
  # Get variance-covariance matrix of reproducible component.
  Sigma <- Reduce(f = "+", x = lapply(seq_along(K), function(i) {
    K[i] * tcrossprod(z[i, ] - mu)
  })) / sum(K)

  # Get variance and standard deviation of reproducible component.
  sigma2 <- mean(diag(Sigma))
  sigma <- sqrt(sigma2)
  
  # Get correlation of reproducible component.
  rho <- mean(Sigma[upper.tri(Sigma)]) / sigma2 
  
  return(list(p = p, 
              mu = mu, 
              sigma = sigma, 
              rho = rho))
}

# Get log likelihood of pseudo data ("Z").
loglihood <- function (z, mu, sigma, rho, p)
{
  l.m <- sum(norm_dens(z, 0, 1, 0) +
               log(p * exp(norm_dens(z, mu, sigma, rho) -
                             norm_dens(z, 0, 1, 0)) +
                     (1 - p)))
  return(l.m)
}

# Define a function to evaluate convergence.
conv <- function(old, new, eps) abs(new - old) < eps * (1 + abs(new))

# Evaluate if the EM algorithm got trapped into an oscillation.
loglik_osc_test <- function(x)
{
  if (length(x) > 10 && 
      stats::var(x[seq(5, length(x), by = 2)]) == 0) stop('Caught in loop!')
}

# Sample eCV values from the null distribution.
ecv_null_sampler <- function(i, n, mu, Sigma, ecv) 
{
  z <- mvtnorm::rmvnorm(n, mean = mu, sigma = Sigma)
  m <- rowMeans(z)
  row_sd <- apply(X = z, MARGIN = 1, FUN = stats::sd, na.rm = TRUE)
  ecv < abs(row_sd^2 - m^2) / m
}
