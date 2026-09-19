\name{debiased_kde}
\alias{debiased_kde}
\alias{kde_confidence_band}
\title{Debiased kernel density estimation and confidence bands}
\description{
Evaluate the Gaussian debiased KDE in equation (3) of Cheng and Chen or
construct its empirical-bootstrap simultaneous confidence band (Figure 2).
}
\usage{
debiased_kde(x, points = NULL, bandwidth = NULL, tau = 1,
  grid_size = 200L)

kde_confidence_band(x, points = NULL, bandwidth = NULL, tau = 1,
  confidence = 0.95, n_boot = 999L, studentized = FALSE,
  random_state = NULL, grid_size = 200L)
}
\arguments{
  \item{x}{Numeric observations or a numeric matrix with observations in rows.}
  \item{points}{Optional evaluation vector or matrix; required for multivariate data.}
  \item{bandwidth}{Positive isotropic bandwidth selected for the ordinary KDE.}
  \item{tau}{Positive ratio \eqn{h/b}; the paper recommends one.}
  \item{grid_size}{Generated grid size for one-dimensional data.}
  \item{confidence}{Confidence level strictly between zero and one.}
  \item{n_boot}{Positive number of empirical-bootstrap replicates.}
  \item{studentized}{Whether to use the variable-width band in Remark 1.}
  \item{random_state}{Optional local random seed; the caller's RNG state is preserved.}
}
\value{
\code{debiased_kde()} returns a \code{di_estimate};
\code{kde_confidence_band()} returns a \code{di_band}.
}
\references{
Cheng, G. and Chen, Y.-C. (2019). Nonparametric Inference via Bootstrapping
the Debiased Estimator. \emph{Electronic Journal of Statistics}, 13(1).
\doi{10.1214/19-EJS1575}.
}
\examples{
x <- rnorm(50)
fit <- debiased_kde(x, bandwidth = 0.4)
band <- kde_confidence_band(x, bandwidth = 0.4, n_boot = 19,
                            random_state = 1)
}
