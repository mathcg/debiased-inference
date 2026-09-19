.crossing_roots <- function(points, values, level) {
  ordering <- order(points)
  x <- points[ordering]
  shifted <- values[ordering] - level
  roots <- x[shifted == 0]
  crossing <- which(
    shifted[seq_len(length(shifted) - 1L)] * shifted[-1L] < 0
  )
  for (index in crossing) {
    fraction <- -shifted[index] / (shifted[index + 1L] - shifted[index])
    roots <- c(roots, x[index] + fraction * (x[index + 1L] - x[index]))
  }
  sort(unique(as.double(roots)))
}

.contour_points <- function(points, values, level) {
  x_coordinates <- sort(unique(points[, 1L]))
  y_coordinates <- sort(unique(points[, 2L]))
  if (nrow(points) != length(x_coordinates) * length(y_coordinates)) {
    stop("two-dimensional level sets require a complete rectangular grid",
         call. = FALSE)
  }
  surface <- matrix(NA_real_, nrow = length(x_coordinates),
                    ncol = length(y_coordinates))
  x_index <- match(points[, 1L], x_coordinates)
  y_index <- match(points[, 2L], y_coordinates)
  surface[cbind(x_index, y_index)] <- values
  if (any(!is.finite(surface))) {
    stop("two-dimensional level sets require unique rectangular grid points",
         call. = FALSE)
  }
  lines <- grDevices::contourLines(
    x = x_coordinates, y = y_coordinates, z = surface, levels = level
  )
  if (length(lines) == 0L) return(matrix(numeric(), nrow = 0L, ncol = 2L))
  geometry <- do.call(rbind, lapply(lines, function(line) cbind(line$x, line$y)))
  unique(round(geometry, digits = 14L))
}

.level_geometry <- function(points, values, level) {
  if (is.null(dim(points))) return(.crossing_roots(points, values, level))
  if (is.matrix(points) && ncol(points) == 2L) {
    return(.contour_points(points, values, level))
  }
  stop("level-set geometry is supported for one or two dimensions", call. = FALSE)
}

.as_set_points <- function(x, name) {
  if (is.vector(x) && is.numeric(x)) x <- matrix(x, ncol = 1L)
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) == 0L || ncol(x) == 0L) {
    stop(name, " must be a non-empty vector or point matrix", call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) stop(name, " must contain only finite values", call. = FALSE)
  x
}

.distance_to_geometry <- function(points, geometry) {
  point_cloud <- if (is.null(dim(points))) matrix(points, ncol = 1L) else points
  geometry_cloud <- if (is.null(dim(geometry))) matrix(geometry, ncol = 1L) else geometry
  apply(point_cloud, 1L, function(point) {
    min(sqrt(rowSums(sweep(geometry_cloud, 2L, point, "-")^2)))
  })
}

#' Estimate a one- or two-dimensional equality level set
#' @param estimate An object returned by `debiased_kde()` or
#'   `debiased_local_linear()`.
#' @param level Finite target level.
#' @return An object of class `di_set`.
#' @export
level_set <- function(estimate, level) {
  if (!inherits(estimate, "di_estimate")) {
    stop("estimate must be a di_estimate object", call. = FALSE)
  }
  points <- estimate$points
  if (is.null(dim(points))) {
    points <- .as_vector(points, "estimate$points")
  } else {
    points <- .as_set_points(points, "estimate$points")
  }
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level)) {
    stop("level must be finite", call. = FALSE)
  }
  structure(list(
    points = points, mask = rep(FALSE, length(estimate$estimate)),
    roots = .level_geometry(points, estimate$estimate, level),
    level = as.double(level), method = "level_set", radius = NULL
  ), class = "di_set")
}

#' Invert a simultaneous confidence band
#'
#' Returns grid points whose intervals contain the target level. This is the
#' alternative density level-set construction in Remark 2 and its inverse-
#' regression analogue.
#' @param band An object returned by a confidence-band function.
#' @param level Finite target level.
#' @return An object of class `di_set`.
#' @export
invert_confidence_band <- function(band, level) {
  if (!inherits(band, "di_band")) stop("band must be a di_band object", call. = FALSE)
  points <- band$points
  if (is.null(dim(points))) {
    points <- .as_vector(points, "band$points")
  } else {
    points <- .as_set_points(points, "band$points")
  }
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level)) {
    stop("level must be finite", call. = FALSE)
  }
  structure(list(
    points = points, mask = band$lower <= level & level <= band$upper,
    roots = .level_geometry(points, band$estimate, level),
    level = as.double(level), method = "band_inversion", radius = NULL,
    confidence = band$confidence, n_boot = band$n_boot,
    bootstrap_statistics = band$bootstrap_statistics
  ), class = "di_set")
}

#' Hausdorff distance between finite point clouds
#' @param a First non-empty numeric vector or point matrix.
#' @param b Second non-empty numeric vector or point matrix.
#' @return A nonnegative scalar.
#' @export
hausdorff_distance <- function(a, b) {
  a <- .as_set_points(a, "a")
  b <- .as_set_points(b, "b")
  if (ncol(a) != ncol(b)) stop("a and b must have the same point dimension", call. = FALSE)
  distances <- matrix(0, nrow = nrow(a), ncol = nrow(b))
  for (coordinate in seq_len(ncol(a))) {
    differences <- outer(a[, coordinate], b[, coordinate], "-")
    distances <- distances + differences^2
  }
  distances <- sqrt(distances)
  max(max(apply(distances, 1L, min)), max(apply(distances, 2L, min)))
}

#' Estimate a density equality level set
#' @inheritParams debiased_kde
#' @param level Finite target density level.
#' @return An object of class `di_set`.
#' @export
density_level_set <- function(x, level, points = NULL, bandwidth = NULL,
                              tau = 1, grid_size = 200L) {
  level_set(debiased_kde(
    x, points, bandwidth = bandwidth, tau = tau, grid_size = grid_size
  ), level)
}

#' Confidence set for a density equality level set
#'
#' The Hausdorff method implements Section 3.1.1 of Cheng and Chen. The
#' inversion method implements the alternative construction in Remark 2.
#' @inheritParams kde_confidence_band
#' @param level Finite target density level.
#' @param method Either `"hausdorff"` or `"inversion"`.
#' @param max_attempts Maximum bootstrap samples attempted if an estimated
#'   level set is empty.
#' @return An object of class `di_set`.
#' @export
density_level_set_confidence <- function(
    x, level, points = NULL, bandwidth = NULL, tau = 1, confidence = 0.95,
    n_boot = 999L, method = "hausdorff", random_state = NULL,
    grid_size = 200L, max_attempts = NULL) {
  if (length(method) != 1L || is.na(method) ||
      !method %in% c("hausdorff", "inversion")) {
    stop("method must be 'hausdorff' or 'inversion'", call. = FALSE)
  }
  if (method == "inversion") {
    return(invert_confidence_band(kde_confidence_band(
      x, points, bandwidth = bandwidth, tau = tau, confidence = confidence,
      n_boot = n_boot, random_state = random_state, grid_size = grid_size
    ), level))
  }
  parameters <- .bootstrap_parameters(confidence, n_boot)
  confidence <- parameters$confidence
  n_boot <- parameters$n_boot
  samples <- .as_samples(x)
  if (!ncol(samples) %in% c(1L, 2L)) {
    stop("Hausdorff level-set confidence supports one or two dimensions",
         call. = FALSE)
  }
  evaluation_matrix <- .as_points(points, samples, grid_size)
  evaluation <- if (ncol(samples) == 1L) {
    evaluation_matrix[, 1L]
  } else {
    evaluation_matrix
  }
  if (is.null(bandwidth)) bandwidth <- density_bandwidth(samples)
  bandwidth <- .positive_scalar(bandwidth, "bandwidth")
  tau <- .positive_scalar(tau, "tau")
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level)) {
    stop("level must be finite", call. = FALSE)
  }
  level <- as.double(level)
  kernel_matrix <- .kernel_matrix(samples, evaluation_matrix, bandwidth, tau)
  estimate <- .kde_values(kernel_matrix, bandwidth, ncol(samples))
  roots <- .level_geometry(evaluation, estimate, level)
  if (length(roots) == 0L) {
    stop("estimated level set is empty on the evaluation grid", call. = FALSE)
  }
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
      indices <- sample.int(nrow(samples), nrow(samples), replace = TRUE)
      counts <- tabulate(indices, nbins = nrow(samples))
      bootstrap_estimate <- as.double(kernel_matrix %*% counts) /
        (nrow(samples) * bandwidth^ncol(samples))
      bootstrap_roots <- .level_geometry(evaluation, bootstrap_estimate, level)
      if (length(bootstrap_roots) > 0L) {
        accepted <- accepted + 1L
        result[accepted] <- hausdorff_distance(bootstrap_roots, roots)
      }
    }
    if (accepted < n_boot) {
      stop("only ", accepted, " non-empty bootstrap level sets after ",
           attempts, " attempts", call. = FALSE)
    }
    result
  })
  radius <- as.double(stats::quantile(
    statistics, confidence, names = FALSE, type = 7
  ))
  distance_to_roots <- .distance_to_geometry(evaluation, roots)
  structure(list(
    points = evaluation, mask = distance_to_roots <= radius, roots = roots,
    level = level, method = "hausdorff", radius = radius,
    confidence = confidence, n_boot = n_boot,
    bootstrap_statistics = statistics
  ), class = "di_set")
}

#' Estimate an inverse-regression equality set
#' @inheritParams debiased_local_linear
#' @param level Finite target response level.
#' @return An object of class `di_set`.
#' @export
inverse_regression <- function(x, y, level, points = NULL, bandwidth = NULL,
                               tau = 1, grid_size = 200L, n_folds = 5L,
                               random_state = 0L) {
  level_set(debiased_local_linear(
    x, y, points, bandwidth = bandwidth, tau = tau, grid_size = grid_size,
    n_folds = n_folds, random_state = random_state
  ), level)
}

#' Confidence set for an inverse-regression equality set
#' @inheritParams regression_confidence_band
#' @param level Finite target response level.
#' @param method One of `"hausdorff"`, `"inversion"`, or `"normal"`. The
#'   normal method requires exactly one estimated crossing.
#' @return An object of class `di_set`.
#' @export
inverse_regression_confidence <- function(
    x, y, level, points = NULL, bandwidth = NULL, tau = 1,
    confidence = 0.95, n_boot = 999L, method = "hausdorff",
    random_state = NULL, grid_size = 200L, n_folds = 5L,
    max_attempts = NULL) {
  if (length(method) != 1L || is.na(method) ||
      !method %in% c("hausdorff", "inversion", "normal")) {
    stop("method must be 'hausdorff', 'inversion', or 'normal'", call. = FALSE)
  }
  if (method == "inversion") {
    return(invert_confidence_band(regression_confidence_band(
      x, y, points, bandwidth = bandwidth, tau = tau,
      confidence = confidence, n_boot = n_boot,
      random_state = random_state, grid_size = grid_size,
      n_folds = n_folds, max_attempts = max_attempts
    ), level))
  }
  parameters <- .bootstrap_parameters(confidence, n_boot)
  confidence <- parameters$confidence
  n_boot <- parameters$n_boot
  x <- .as_vector(x, "x")
  y <- .as_vector(y, "y")
  if (length(x) != length(y) || length(x) < 4L) {
    stop("x and y must have equal lengths of at least 4", call. = FALSE)
  }
  base <- debiased_local_linear(
    x, y, points, bandwidth = bandwidth, tau = tau, grid_size = grid_size,
    n_folds = n_folds, random_state = random_state
  )
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level)) {
    stop("level must be finite", call. = FALSE)
  }
  level <- as.double(level)
  roots <- .crossing_roots(base$points, base$estimate, level)
  if (length(roots) == 0L) {
    stop("estimated inverse-regression set is empty on the evaluation grid",
         call. = FALSE)
  }
  if (method == "normal" && length(roots) != 1L) {
    stop("method='normal' requires exactly one estimated crossing", call. = FALSE)
  }
  if (method == "normal" && n_boot < 2L) {
    stop("method='normal' requires at least two bootstrap replicates", call. = FALSE)
  }
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
        bootstrap_roots <- .crossing_roots(base$points, estimate, level)
        valid_roots <- if (method == "normal") {
          length(bootstrap_roots) == 1L
        } else {
          length(bootstrap_roots) > 0L
        }
        if (valid_roots) {
          accepted <- accepted + 1L
          result[accepted] <- if (method == "normal") {
            bootstrap_roots[1L]
          } else {
            hausdorff_distance(bootstrap_roots, roots)
          }
        }
      }
    }
    if (accepted < n_boot) {
      stop("only ", accepted, " non-empty bootstrap inverse sets after ",
           attempts, " attempts", call. = FALSE)
    }
    result
  })
  radius <- if (method == "normal") {
    stats::qnorm(0.5 + confidence / 2) * stats::sd(statistics)
  } else {
    as.double(stats::quantile(statistics, confidence, names = FALSE, type = 7))
  }
  distance_to_roots <- apply(abs(outer(base$points, roots, "-")), 1L, min)
  structure(list(
    points = base$points, mask = distance_to_roots <= radius, roots = roots,
    level = level, method = method, radius = radius,
    confidence = confidence, n_boot = n_boot,
    bootstrap_statistics = statistics
  ), class = "di_set")
}
