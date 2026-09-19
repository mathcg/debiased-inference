\name{debiased_local_linear}
\alias{debiased_local_linear}
\alias{regression_confidence_band}
\title{Debiased local-linear regression and confidence bands}
\description{
Fit the debiased local-linear estimator in equation (5) of Cheng and Chen or
construct its paired-bootstrap simultaneous confidence band (Figure 3).
}
\usage{
debiased_local_linear(x, y, points = NULL, bandwidth = NULL, tau = 1,
  grid_size = 200L, n_folds = 5L, random_state = 0L)

regression_confidence_band(x, y, points = NULL, bandwidth = NULL,
  tau = 1, confidence = 0.95, n_boot = 999L, random_state = NULL,
  grid_size = 200L, n_folds = 5L, max_attempts = NULL)
}
\arguments{
  \item{x}{Numeric one-dimensional covariates.}
  \item{y}{Numeric responses of the same length as \code{x}.}
  \item{points}{Optional numeric evaluation grid.}
  \item{bandwidth}{Positive bandwidth selected for the ordinary smoother.}
  \item{tau}{Positive ratio \eqn{h/b}; the paper recommends one.}
  \item{grid_size}{Generated grid size when \code{points} is omitted.}
  \item{n_folds}{Number of cross-validation folds for bandwidth selection.}
  \item{random_state}{Optional local random seed.}
  \item{confidence}{Confidence level strictly between zero and one.}
  \item{n_boot}{Positive number of paired-bootstrap replicates.}
  \item{max_attempts}{Maximum resamples attempted if singular fits occur.}
}
\value{A \code{di_estimate} or \code{di_band} object.}
\references{
Cheng, G. and Chen, Y.-C. (2019). Nonparametric Inference via Bootstrapping
the Debiased Estimator. \emph{Electronic Journal of Statistics}, 13(1).
\doi{10.1214/19-EJS1575}.
}
\examples{
x <- seq(-1, 1, length.out = 50)
y <- sin(pi * x) + rnorm(50, sd = 0.1)
fit <- debiased_local_linear(x, y, bandwidth = 0.3)
band <- regression_confidence_band(x, y, bandwidth = 0.3,
                                   n_boot = 19, random_state = 1)
}
