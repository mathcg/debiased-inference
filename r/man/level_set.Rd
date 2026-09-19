\name{level_set}
\alias{level_set}
\alias{invert_confidence_band}
\alias{hausdorff_distance}
\alias{density_level_set}
\alias{density_level_set_confidence}
\alias{inverse_regression}
\alias{inverse_regression_confidence}
\title{Grid-based level-set utilities}
\description{
Estimate equality level sets on a one-dimensional or rectangular two-
dimensional grid, invert a simultaneous band to obtain a confidence set, or
calculate finite-point-cloud Hausdorff distance.
}
\usage{
level_set(estimate, level)

invert_confidence_band(band, level)

hausdorff_distance(a, b)

density_level_set(x, level, points = NULL, bandwidth = NULL, tau = 1,
  grid_size = 200L)

density_level_set_confidence(x, level, points = NULL, bandwidth = NULL,
  tau = 1, confidence = 0.95, n_boot = 999L, method = "hausdorff",
  random_state = NULL, grid_size = 200L, max_attempts = NULL)

inverse_regression(x, y, level, points = NULL, bandwidth = NULL,
  tau = 1, grid_size = 200L, n_folds = 5L, random_state = 0L)

inverse_regression_confidence(x, y, level, points = NULL,
  bandwidth = NULL, tau = 1, confidence = 0.95, n_boot = 999L,
  method = "hausdorff", random_state = NULL, grid_size = 200L,
  n_folds = 5L, max_attempts = NULL)
}
\arguments{
  \item{estimate}{A \code{di_estimate} object.}
  \item{band}{A \code{di_band} object.}
  \item{level}{Finite target level.}
  \item{a}{First non-empty numeric vector or point matrix.}
  \item{b}{Second non-empty numeric vector or point matrix.}
  \item{x}{Numeric observations or covariates.}
  \item{y}{Numeric responses for inverse regression.}
  \item{points}{Optional evaluation grid. Density level sets accept a vector
  or a complete rectangular two-dimensional point matrix; inverse regression
  accepts a vector.}
  \item{bandwidth}{Optional positive bandwidth.}
  \item{tau}{Positive ratio \eqn{h/b}.}
  \item{grid_size}{Generated grid size when points are omitted.}
  \item{confidence}{Confidence level strictly between zero and one.}
  \item{n_boot}{Positive number of bootstrap replicates.}
  \item{method}{Either \code{"hausdorff"} or \code{"inversion"}; inverse
  regression additionally supports \code{"normal"} for a unique crossing.}
  \item{random_state}{Optional local random seed.}
  \item{max_attempts}{Maximum resamples attempted when sets are empty.}
  \item{n_folds}{Cross-validation folds for a regression bandwidth.}
}
\value{Level-set and inverse-regression functions return a \code{di_set};
\code{hausdorff_distance()} returns a nonnegative scalar.}
