.local_polynomial <- function(x, y, points, bandwidth, degree, derivative) {
  result <- rep(NA_real_, length(points))
  for (point_index in seq_along(points)) {
    offsets <- x - points[point_index]
    weights <- stats::dnorm(offsets / bandwidth)
    design <- outer(offsets, 0:degree, "^")
    root_weights <- sqrt(weights)
    weighted_design <- design * root_weights
    weighted_response <- y * root_weights
    fit <- qr(weighted_design, tol = 1e-10)
    if (fit$rank == degree + 1L) {
      coefficients <- qr.coef(fit, weighted_response)
      result[point_index] <- factorial(derivative) * coefficients[derivative + 1L]
    }
  }
  result
}

.debiased_regression_values <- function(x, y, points, bandwidth, tau) {
  ordinary <- .local_polynomial(
    x, y, points, bandwidth, degree = 1L, derivative = 0L
  )
  second <- .local_polynomial(
    x, y, points, bandwidth / tau, degree = 3L, derivative = 2L
  )
  ordinary - 0.5 * bandwidth^2 * second
}

#' Debiased local-linear regression estimator
#'
#' Fits the ordinary local-linear smoother and subtracts one half h-squared
#' times a local-cubic estimate of the second derivative, as in equation (5)
#' of Cheng and Chen.
#' @param x Numeric one-dimensional covariates.
#' @param y Numeric responses.
#' @param points Optional numeric evaluation grid.
#' @param bandwidth Positive scalar selected for the ordinary smoother.
#' @param tau Positive ratio h/b; the paper recommends 1.
#' @param grid_size Generated grid size when `points` is omitted.
#' @param n_folds Folds used when selecting a bandwidth.
#' @param random_state Optional local random seed.
#' @return An object of class `di_estimate`.
#' @export
debiased_local_linear <- function(x, y, points = NULL, bandwidth = NULL,
                                  tau = 1, grid_size = 200L, n_folds = 5L,
                                  random_state = 0L) {
  x <- .as_vector(x, "x")
  y <- .as_vector(y, "y")
  if (length(x) != length(y) || length(x) < 4L) {
    stop("x and y must have equal lengths of at least 4", call. = FALSE)
  }
  if (diff(range(x)) == 0) stop("x must contain at least two distinct values", call. = FALSE)
  if (is.null(points)) {
    if (!is.numeric(grid_size) || length(grid_size) != 1L || grid_size < 2 ||
        grid_size != as.integer(grid_size)) {
      stop("grid_size must be an integer of at least 2", call. = FALSE)
    }
    points <- seq(min(x), max(x), length.out = as.integer(grid_size))
  } else {
    points <- .as_vector(points, "points")
  }
  if (is.null(bandwidth)) {
    bandwidth <- regression_bandwidth(
      x, y, n_folds = n_folds, random_state = random_state
    )
  }
  bandwidth <- .positive_scalar(bandwidth, "bandwidth")
  tau <- .positive_scalar(tau, "tau")
  estimate <- .debiased_regression_values(x, y, points, bandwidth, tau)
  if (any(!is.finite(estimate))) {
    stop("singular local fit; increase bandwidth or restrict evaluation points",
         call. = FALSE)
  }
  structure(list(
    points = points, estimate = estimate, bandwidth = bandwidth, tau = tau,
    method = "debiased_local_linear"
  ), class = "di_estimate")
}

#' Simultaneous paired-bootstrap confidence band for a regression function
#'
#' Implements Figure 3 of Cheng and Chen.
#' @inheritParams debiased_local_linear
#' @param confidence Confidence level strictly between zero and one.
#' @param n_boot Number of paired-bootstrap replicates.
#' @param max_attempts Maximum resamples attempted when singular fits are rejected.
#' @return An object of class `di_band`.
#' @export
regression_confidence_band <- function(
    x, y, points = NULL, bandwidth = NULL, tau = 1, confidence = 0.95,
    n_boot = 999L, random_state = NULL, grid_size = 200L, n_folds = 5L,
    max_attempts = NULL) {
  parameters <- .bootstrap_parameters(confidence, n_boot)
  confidence <- parameters$confidence
  n_boot <- parameters$n_boot
  base <- debiased_local_linear(
    x, y, points, bandwidth = bandwidth, tau = tau, grid_size = grid_size,
    n_folds = n_folds, random_state = random_state
  )
  x <- .as_vector(x, "x")
  y <- .as_vector(y, "y")
  if (is.null(max_attempts)) max_attempts <- max(10L * n_boot, 100L)
  if (!is.numeric(max_attempts) || length(max_attempts) != 1L ||
      !is.finite(max_attempts) ||
      max_attempts < n_boot || max_attempts != as.integer(max_attempts)) {
    stop("max_attempts must be an integer at least n_boot", call. = FALSE)
  }
  statistics <- .with_seed(random_state, {
    result <- numeric(n_boot)
    accepted <- 0L
    attempts <- 0L
    while (accepted < n_boot && attempts < max_attempts) {
      attempts <- attempts + 1L
      indices <- sample.int(length(x), length(x), replace = TRUE)
      estimate <- .debiased_regression_values(
        x[indices], y[indices], base$points, base$bandwidth, base$tau
      )
      if (all(is.finite(estimate))) {
        accepted <- accepted + 1L
        result[accepted] <- max(abs(estimate - base$estimate))
      }
    }
    if (accepted < n_boot) {
      stop("only ", accepted, " nonsingular bootstrap fits after ", attempts,
           " attempts", call. = FALSE)
    }
    result
  })
  critical_value <- as.double(stats::quantile(
    statistics, confidence, names = FALSE, type = 7
  ))
  structure(list(
    points = base$points, estimate = base$estimate,
    lower = base$estimate - critical_value,
    upper = base$estimate + critical_value,
    critical_value = critical_value, bandwidth = base$bandwidth, tau = base$tau,
    confidence = confidence, n_boot = n_boot,
    bootstrap_statistics = statistics, method = "debiased_local_linear",
    studentized = FALSE, standard_error = NULL
  ), class = "di_band")
}
