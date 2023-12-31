% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gIDR.R
\name{gIDR}
\alias{gIDR}
\title{Estimate IDR with a General Mixture Model}
\usage{
gIDR(x, mu, sigma, rho, p, eps = 0.001, max.ite = 30)
}
\arguments{
\item{x}{Numeric matrix with rows representing the number of omic features 
and columns representing the number of sample replicates.
The numeric values must be positive and represent significance 
(not necessarily p-values).}

\item{mu}{Starting value for the mean of the reproducible component. Numeric.}

\item{sigma}{Starting value for the standard deviation of the reproducible 
component.}

\item{rho}{Starting value for the correlation coefficient of the reproducible 
component.}

\item{p}{Starting value for the proportion of the reproducible component.}

\item{eps}{Stopping criterion. Iterations stop when the increment of the 
log-likelihood is less than eps times the log-likelihood. Default is 0.001.}

\item{max.ite}{Maximum number of iterations. Default is 30.}
}
\value{
Returns a list of three elements:
\describe{
  \item{\strong{idr}}{A numeric vector of the local IDR 
  (Irreproducible Discovery Rate) for each observation (i.e., estimated 
  conditional probability for each observation to belong to the 
  irreproducible component).}
  \item{\strong{IDR}}{A numerical vector of the expected Irreproducible
   Discovery Rate for observations that are as irreproducible or more 
   irreproducible than the given observations.}
  \item{\strong{est_param}}{Estimated parameters: p, rho, mu, sigma.}

}
}
\description{
This function builds  upon `idr::est.IDR` to extend the Li et al. (2011)
copula mixture model to accommodate  an  arbitrary number of replicates.
The  term  "General"  in  this context  alludes to  the assumption of a 
general multivariate Normal distribution of dimension "n," equal to the
number  of  sample  replicates. This  assumption essentially allows the
pseudo-likelihood  approach  in  Li et al. (2011) to be extended to any
number of replicates.  This  is  achieved  by  modifying the "E" and "M" 
steps  of  an  expectation  maximization algorithm to use a multivariate
Normal Distribution instead.
}
\examples{
# 1. Show that gIDR reduces to classical IDR for n=2.

# Load required packages.
library(idr)
library(eCV)

# Set seed for RNG.
set.seed(42)

# Simulate data.
out <- simulate_data(scenario = 2, n_features = 1e3)

# Set initial parameter values.
mu <- 2
sigma <- 1.3
rho <- 0.8
p <- 0.7

# Compare IDR and gIDR
idr.out <- est.IDR(x = out$sim_data, mu, sigma, rho, p)
gidr.out <- gIDR(x = out$sim_data, mu, sigma, rho, p)

# Show the results are the same.
all.equal(gidr.out$est_param, idr.out$para)

\donttest{

library(tidyverse)
# Plot results.
out$sim_data \%>\% as.data.frame() \%>\% 
mutate(idr = gidr.out$idr) \%>\%
ggplot(aes(x=`Rep 1`,y=`Rep 2`,color=idr)) +
 geom_point(size=1) + 
 scale_color_gradientn(colors=c("#F4364C", "#D5DADD", "#009CA6" ))+ 
 theme_classic()

#2. Show gIDR for n=10.

out <- simulate_data(scenario = 1, n_reps = 10, n_features = 1e3)
gidr.out <- gIDR(x = out$sim_data, mu, sigma, rho, p)
out$sim_data \%>\% as.data.frame() \%>\% 
mutate(idr = gidr.out$IDR) \%>\%
ggplot(aes(x = `Rep 1`, y = `Rep 2`, color = idr)) +
 geom_point(size = 1) + 
 scale_color_gradientn(colors=c("#F4364C", "#D5DADD", "#009CA6" ))+ 
 theme_classic()
}
}
\references{
Q. Li, J. B. Brown, H. Huang, and P. J. Bickel. (2011)
}
