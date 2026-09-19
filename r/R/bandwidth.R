#' Select an ordinary density-estimation bandwidth
#'
#' Uses the Gaussian normal-scale rule used by the paper's `ks`
#' implementation. In one dimension this is `(4 / (3 n))^(1/5) * sd`.
#' Bandwidth selection is for the ordinary estimator, as required by the paper.
#' @param x Numeric observations or a numeric matrix with observations in rows.
#' @param method Either `"normal_reference"` or `"cv"`.
#' @param candidates Optional positive candidate bandwidths for cross-validation.
#' @param block_size Pairwise-computation block size for density cross-validation.
#' @return A positive numeric scalar.
#' @export
density_bandwidth <- function(x, method = "normal_reference", candidates = NULL,
                              block_size = 512L) {
  samples <- .as_samples(x)
  if (length(method) != 1L || is.na(method) ||
      !method %in% c("normal_reference", "cv")) {
    stop("method must be 'normal_reference' or 'cv'", call. = FALSE)
  }
  n <- nrow(samples)
  dimension <- ncol(samples)
  standard_deviation <- apply(samples, 2L, stats::sd)
  positive <- standard_deviation[standard_deviation > 0]
  if (length(positive) == 0L) {
    stop("cannot select bandwidth from zero-scale data", call. = FALSE)
  }
  scale <- exp(mean(log(positive)))
  factor <- (4 / (dimension + 2))^(1 / (dimension + 4)) *
    n^(-1 / (dimension + 4))
  reference <- .positive_scalar(scale * factor, "selected bandwidth")
  if (method == "normal_reference") {
    if (!is.null(candidates)) {
      stop("candidates are only used when method='cv'", call. = FALSE)
    }
    return(reference)
  }
  if (is.null(candidates)) {
    candidates <- reference * exp(seq(log(0.35), log(2.5), length.out = 31L))
  } else {
    candidates <- .as_vector(candidates, "candidates")
    if (any(candidates <= 0)) stop("candidates must be positive", call. = FALSE)
  }
  if (!is.numeric(block_size) || length(block_size) != 1L ||
      !is.finite(block_size) || block_size < 1 ||
      block_size != as.integer(block_size)) {
    stop("block_size must be a positive integer", call. = FALSE)
  }
  block_size <- as.integer(block_size)
  scores <- numeric(length(candidates))
  for (index in seq_along(candidates)) {
    h <- candidates[index]
    integrated_normalizer <- (4 * pi * h^2)^(-dimension / 2)
    ordinary_normalizer <- (2 * pi * h^2)^(-dimension / 2)
    integrated_sum <- 0
    ordinary_sum <- 0
    starts <- seq.int(1L, n, by = block_size)
    for (start in starts) {
      stop_index <- min(start + block_size - 1L, n)
      block <- samples[start:stop_index, , drop = FALSE]
      norm_squared <- matrix(0, nrow = nrow(block), ncol = n)
      for (coordinate in seq_len(dimension)) {
        differences <- outer(block[, coordinate], samples[, coordinate], "-")
        norm_squared <- norm_squared + differences^2
      }
      integrated_sum <- integrated_sum + integrated_normalizer *
        sum(exp(-norm_squared / (4 * h^2)))
      ordinary_sum <- ordinary_sum + ordinary_normalizer *
        sum(exp(-norm_squared / (2 * h^2)))
    }
    ordinary_off_diagonal <- ordinary_sum - n * ordinary_normalizer
    scores[index] <- integrated_sum / n^2 -
      2 * ordinary_off_diagonal / (n * (n - 1))
  }
  as.double(candidates[which.min(scores)])
}

#' Select a local-linear regression bandwidth
#'
#' Performs K-fold cross-validation on the ordinary local-linear estimator.
#' @param x Numeric covariates.
#' @param y Numeric responses.
#' @param candidates Optional positive candidate bandwidths.
#' @param n_folds Number of folds.
#' @param random_state Optional local random seed.
#' @return The candidate with smallest mean squared validation error.
#' @export
regression_bandwidth <- function(x, y, candidates = NULL, n_folds = 5L,
                                 random_state = 0L) {
  x <- .as_vector(x, "x")
  y <- .as_vector(y, "y")
  if (length(x) != length(y)) stop("x and y must have the same length", call. = FALSE)
  if (diff(range(x)) == 0) stop("x must contain at least two distinct values", call. = FALSE)
  n <- length(x)
  if (!is.numeric(n_folds) || length(n_folds) != 1L ||
      n_folds < 2 || n_folds > n || n_folds != as.integer(n_folds)) {
    stop("n_folds must be an integer between 2 and length(x)", call. = FALSE)
  }
  n_folds <- as.integer(n_folds)
  if (is.null(candidates)) {
    scale <- min(stats::sd(x), diff(range(x)) / 4)
    base <- max(scale * n^(-1 / 5), .Machine$double.eps)
    # Include the small CV bandwidths selected for the paper's high-curvature
    # regression designs; a narrower range can truncate the optimum.
    candidates <- base * exp(seq(log(0.1), log(3), length.out = 31L))
  } else {
    candidates <- .as_vector(candidates, "candidates")
    if (any(candidates <= 0)) stop("candidates must be positive", call. = FALSE)
  }
  fold_id <- .with_seed(random_state, {
    order <- sample.int(n)
    result <- integer(n)
    result[order] <- rep(seq_len(n_folds), length.out = n)
    result
  })
  losses <- rep(Inf, length(candidates))
  for (candidate_index in seq_along(candidates)) {
    errors <- numeric()
    valid <- TRUE
    for (fold in seq_len(n_folds)) {
      held_out <- fold_id == fold
      predicted <- .local_polynomial(
        x[!held_out], y[!held_out], x[held_out], candidates[candidate_index],
        degree = 1L, derivative = 0L
      )
      if (any(!is.finite(predicted))) {
        valid <- FALSE
        break
      }
      errors <- c(errors, (y[held_out] - predicted)^2)
    }
    if (valid) losses[candidate_index] <- mean(errors)
  }
  if (!any(is.finite(losses))) {
    stop("all candidate bandwidths produced singular local fits", call. = FALSE)
  }
  as.double(candidates[which.min(losses)])
}
