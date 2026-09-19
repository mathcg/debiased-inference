\name{density_bandwidth}
\alias{density_bandwidth}
\alias{regression_bandwidth}
\title{Bandwidth selectors for ordinary nonparametric estimators}
\description{
Select the bandwidth on the ordinary estimator, as required by the debiasing
method. Density estimation uses a normal-reference rule and regression uses
K-fold cross-validation of the ordinary local-linear smoother.
}
\usage{
density_bandwidth(x, method = "normal_reference", candidates = NULL,
  block_size = 512L)

regression_bandwidth(x, y, candidates = NULL, n_folds = 5L,
  random_state = 0L)
}
\arguments{
  \item{x}{Numeric observations or covariates.}
  \item{y}{Numeric responses.}
  \item{method}{For density estimation, either \code{"normal_reference"} or
  \code{"cv"}.}
  \item{candidates}{Optional positive candidate bandwidths for cross-validation.}
  \item{block_size}{Positive pairwise-computation block size for density
  cross-validation.}
  \item{n_folds}{Number of cross-validation folds.}
  \item{random_state}{Optional local random seed.}
}
\value{A positive numeric scalar.}
