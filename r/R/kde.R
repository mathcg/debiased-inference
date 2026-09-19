.gaussian_debiased_kernel <- function(norm_squared, dimension, tau) {
  normalizer <- (2 * pi)^(-dimension / 2)
  ordinary <- normalizer * exp(-0.5 * norm_squared)
  tau_norm_squared <- tau^2 * norm_squared
  laplacian <- normalizer * exp(-0.5 * tau_norm_squared) *
    (tau_norm_squared - dimension)
  ordinary - 0.5 * tau^(dimension + 2) * laplacian
}

.kernel_matrix <- function(samples, points, bandwidth, tau) {
  norm_squared <- matrix(0, nrow = nrow(points), ncol = nrow(samples))
  for (coordinate in seq_len(ncol(samples))) {
    differences <- outer(points[, coordinate], samples[, coordinate], "-") / bandwidth
    norm_squared <- norm_squared + differences^2
  }
  .gaussian_debiased_kernel(norm_squared, ncol(samples), tau)
}

.kde_values <- function(kernel_matrix, bandwidth, dimension) {
  rowMeans(kernel_matrix) / bandwidth^dimension
}

#' Debiased kernel density estimator
#'
#' Evaluates equation (3) of Cheng and Chen using the Gaussian kernel.
#' @param x Numeric observations or a numeric matrix with observations in rows.
#' @param points Optional evaluation vector or matrix. Required for multivariate data.
#' @param bandwidth Positive scalar selected for the ordinary KDE.
#' @param tau Positive ratio h/b; the paper recommends 1.
#' @param grid_size Generated grid size for one-dimensional data.
#' @return An object of class `di_estimate`.
#' @export
debiased_kde <- function(x, points = NULL, bandwidth = NULL, tau = 1,
                         grid_size = 200L) {
  samples <- .as_samples(x)
  evaluation <- .as_points(points, samples, grid_size)
  if (is.null(bandwidth)) bandwidth <- density_bandwidth(samples)
  bandwidth <- .positive_scalar(bandwidth, "bandwidth")
  tau <- .positive_scalar(tau, "tau")
  kernel_matrix <- .kernel_matrix(samples, evaluation, bandwidth, tau)
  structure(list(
    points = .public_points(evaluation),
    estimate = .kde_values(kernel_matrix, bandwidth, ncol(samples)),
    bandwidth = bandwidth,
    tau = tau,
    method = "debiased_kde"
  ), class = "di_estimate")
}

#' Simultaneous bootstrap confidence band for a density
#'
#' Implements Figure 2 of Cheng and Chen. Setting `studentized = TRUE`
#' implements the variable-width band in Remark 1.
#' @inheritParams debiased_kde
#' @param confidence Confidence level strictly between zero and one.
#' @param n_boot Number of empirical-bootstrap replicates.
#' @param studentized Whether to construct the studentized variable-width band.
#' @param random_state Optional local random seed.
#' @return An object of class `di_band`.
#' @export
kde_confidence_band <- function(x, points = NULL, bandwidth = NULL, tau = 1,
                                confidence = 0.95, n_boot = 999L,
                                studentized = FALSE, random_state = NULL,
                                grid_size = 200L) {
  parameters <- .bootstrap_parameters(confidence, n_boot)
  confidence <- parameters$confidence
  n_boot <- parameters$n_boot
  if (!is.logical(studentized) || length(studentized) != 1L || is.na(studentized)) {
    stop("studentized must be TRUE or FALSE", call. = FALSE)
  }
  samples <- .as_samples(x)
  evaluation <- .as_points(points, samples, grid_size)
  if (is.null(bandwidth)) bandwidth <- density_bandwidth(samples)
  bandwidth <- .positive_scalar(bandwidth, "bandwidth")
  tau <- .positive_scalar(tau, "tau")
  dimension <- ncol(samples)
  kernel_matrix <- .kernel_matrix(samples, evaluation, bandwidth, tau)
  estimate <- .kde_values(kernel_matrix, bandwidth, dimension)

  squared_kernel_matrix <- kernel_matrix^2
  first_moment <- rowMeans(kernel_matrix)
  second_moment <- rowMeans(squared_kernel_matrix)
  sigma_squared <- pmax((second_moment - first_moment^2) / bandwidth^dimension, 0)
  standard_error <- sqrt(sigma_squared / (nrow(samples) * bandwidth^dimension))
  if (studentized && any(standard_error <= .Machine$double.eps)) {
    stop("studentized band is undefined where estimated variance is zero", call. = FALSE)
  }
  statistics <- .with_seed(random_state, {
    result <- numeric(n_boot)
    for (bootstrap_index in seq_len(n_boot)) {
      indices <- sample.int(nrow(samples), nrow(samples), replace = TRUE)
      counts <- tabulate(indices, nbins = nrow(samples))
      bootstrap_estimate <- as.double(kernel_matrix %*% counts) /
        (nrow(samples) * bandwidth^dimension)
      difference <- abs(bootstrap_estimate - estimate)
      if (studentized) {
        boot_first <- as.double(kernel_matrix %*% counts) / nrow(samples)
        boot_second <- as.double(squared_kernel_matrix %*% counts) / nrow(samples)
        boot_sigma_squared <- pmax(
          (boot_second - boot_first^2) / bandwidth^dimension, 0
        )
        bootstrap_se <- sqrt(
          boot_sigma_squared / (nrow(samples) * bandwidth^dimension)
        )
        if (any(bootstrap_se <= .Machine$double.eps)) {
          stop("a bootstrap resample has zero estimated variance", call. = FALSE)
        }
        difference <- difference / bootstrap_se
      }
      result[bootstrap_index] <- max(difference)
    }
    result
  })
  critical_value <- as.double(stats::quantile(
    statistics, confidence, names = FALSE, type = 7
  ))
  half_width <- if (studentized) {
    critical_value * standard_error
  } else {
    rep(critical_value, length(estimate))
  }
  structure(list(
    points = .public_points(evaluation), estimate = estimate,
    lower = estimate - half_width, upper = estimate + half_width,
    critical_value = critical_value, bandwidth = bandwidth, tau = tau,
    confidence = confidence, n_boot = n_boot,
    bootstrap_statistics = statistics, method = "debiased_kde",
    studentized = studentized,
    standard_error = if (studentized) standard_error else NULL
  ), class = "di_band")
}
