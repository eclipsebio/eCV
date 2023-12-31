% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eCV.R
\name{eCV}
\alias{eCV}
\title{Assess Reproducibility via Enhanced Coefficient of Variation}
\usage{
eCV(x, max.ite = 10000, n_threads = 1)
}
\arguments{
\item{x}{A numeric matrix with rows representing omic features and columns 
representing sample replicates. Numeric values should be positive and reflect 
significance (not necessarily p-values).}

\item{max.ite}{Number of samples from the null distribution (numeric). 
Defaults to 1e4.}

\item{n_threads}{Number of threads for parallel computing. Numeric. Defaults 
to 1.}
}
\value{
Returns a list with two elements:
\itemize{
\item \emph{\strong{ecv:}} Numeric vector with
the estimated eCV values for each omic feature.
\item \emph{\strong{post_prob:}} Posterior probability values.
}
}
\description{
This function estimates an "enhanced" coefficient of variation (eCV) to 
measure the likelihood of an omic feature being reproducible.  The  eCV 
is calculated as |SD^2 - Mean^2| / Mean^2, a metric that decreases with
noise among replicates and increases with the mean intensity.
}
\details{
Inferences are made based on the probabilities of eCV values originating
from the group of reproducible features.  It assumes  that  reproducible
features follow a prior Normal distribution with dimension r (number  of
replicates). Pseudo replicates are sampled using a Probabilistic Bootstrap,
assuming that the global mean vector and variance-covariance matrix across
features are close to the prior's hyperparameters. 
eCV is computed for each random sample. 
The proportion of times the observed eCV is lower than or equal to the eCV
from random samples is then taken as the probability of  the  omic feature
belonging to the group of reproducible features.
}
\examples{
library(eCV)
set.seed(42)
# Simulate data.
out <- simulate_data(scenario = 1,n_features=1e3)

# Run eCV
ecv_out <- eCV(x = out$sim_data, max.ite = 100)

\donttest{
# Plot results.
library(tidyverse)

out$sim_data \%>\% 
as.data.frame() \%>\% 
mutate(`eCV Prob` = ecv_out$post_prob) \%>\%
ggplot(aes(x = `Rep 1`, y = `Rep 2`, color = `eCV Prob`)) +
 geom_point(size = 1) + 
 scale_color_gradientn(colors=c( "#009CA6", "#D5DADD", "#F4364C"))+ 
 theme_classic()
 }
}
