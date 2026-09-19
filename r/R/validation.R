.as_samples <- function(x, name = "x") {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.vector(x) && is.numeric(x)) x <- matrix(x, ncol = 1L)
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) < 2L || ncol(x) < 1L) {
    stop(name, " must contain at least two observations", call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) {
    stop(name, " must contain only finite values", call. = FALSE)
  }
  ranges <- apply(x, 2L, function(column) diff(range(column)))
  if (all(ranges == 0)) {
    stop(name, " must contain at least two distinct observations", call. = FALSE)
  }
  x
}

.as_vector <- function(x, name) {
  if (!is.numeric(x) || !is.atomic(x) || !is.null(dim(x)) || length(x) == 0L) {
    stop(name, " must be a non-empty numeric vector", call. = FALSE)
  }
  x <- as.double(x)
  if (any(!is.finite(x))) {
    stop(name, " must contain only finite values", call. = FALSE)
  }
  x
}

.as_points <- function(points, samples, grid_size = 200L) {
  dimension <- ncol(samples)
  if (is.null(points)) {
    if (dimension != 1L) {
      stop("points are required for multivariate data", call. = FALSE)
    }
    if (length(grid_size) != 1L || !is.finite(grid_size) ||
        grid_size < 2 || grid_size != as.integer(grid_size)) {
      stop("grid_size must be an integer of at least 2", call. = FALSE)
    }
    spread <- stats::sd(samples[, 1L])
    padding <- 0.05 * max(diff(range(samples[, 1L])), spread)
    if (padding == 0) padding <- 1
    return(matrix(seq(
      min(samples[, 1L]) - padding,
      max(samples[, 1L]) + padding,
      length.out = as.integer(grid_size)
    ), ncol = 1L))
  }
  if (dimension == 1L && is.vector(points) && is.numeric(points)) {
    points <- matrix(points, ncol = 1L)
  }
  if (is.data.frame(points)) points <- as.matrix(points)
  if (!is.matrix(points) || !is.numeric(points) || nrow(points) == 0L ||
      ncol(points) != dimension || any(!is.finite(points))) {
    stop("points must be a finite numeric matrix with ", dimension,
         " column(s)", call. = FALSE)
  }
  storage.mode(points) <- "double"
  points
}

.positive_scalar <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L ||
      !is.finite(value) || value <= 0) {
    stop(name, " must be a finite positive scalar", call. = FALSE)
  }
  as.double(value)
}

.bootstrap_parameters <- function(confidence, n_boot) {
  if (!is.numeric(confidence) || length(confidence) != 1L ||
      !is.finite(confidence) || confidence <= 0 || confidence >= 1) {
    stop("confidence must be strictly between 0 and 1", call. = FALSE)
  }
  if (!is.numeric(n_boot) || length(n_boot) != 1L || !is.finite(n_boot) ||
      n_boot < 1 || n_boot != as.integer(n_boot)) {
    stop("n_boot must be a positive integer", call. = FALSE)
  }
  list(confidence = as.double(confidence), n_boot = as.integer(n_boot))
}

.with_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
    stop("random_state must be NULL or a finite numeric scalar", call. = FALSE)
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(code)
}

.public_points <- function(points) {
  if (ncol(points) == 1L) as.double(points[, 1L]) else points
}
