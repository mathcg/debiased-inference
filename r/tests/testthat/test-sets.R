test_that("level-set crossings are linearly interpolated", {
  estimate <- structure(list(
    points = c(-1, 0, 1), estimate = c(1, -1, 1),
    bandwidth = 0.2, tau = 1, method = "test"
  ), class = "di_estimate")
  result <- level_set(estimate, 0)
  expect_equal(result$roots, c(-0.5, 0.5))
})

test_that("Hausdorff distance is correct", {
  expect_equal(hausdorff_distance(c(0, 2), c(0.5, 2.5)), 0.5)
  expect_error(hausdorff_distance(numeric(), 1))
})

test_that("two-dimensional contours and Hausdorff distance work", {
  axis <- seq(-1, 1, length.out = 9)
  points <- as.matrix(expand.grid(axis, axis))
  estimate <- structure(list(
    points = points,
    estimate = points[, 1]^2 + points[, 2]^2,
    bandwidth = 0.2, tau = 1, method = "test"
  ), class = "di_estimate")
  result <- level_set(estimate, 0.5)
  expect_equal(ncol(result$roots), 2)
  expect_gt(nrow(result$roots), 8)
  expect_equal(hausdorff_distance(result$roots, result$roots), 0)
})

test_that("two-dimensional density level-set confidence works", {
  set.seed(7)
  sample <- matrix(rnorm(120), ncol = 2)
  axis <- seq(-2.5, 2.5, length.out = 17)
  points <- as.matrix(expand.grid(axis, axis))
  result <- density_level_set_confidence(
    sample, 0.05, points, bandwidth = 0.7, n_boot = 5, random_state = 4
  )
  expect_equal(ncol(result$roots), 2)
  expect_length(result$mask, nrow(points))
})

test_that("density level-set Hausdorff confidence is reproducible", {
  x <- c(seq(-1.2, -0.3, length.out = 20), seq(0.3, 1.2, length.out = 20))
  points <- seq(-2, 2, length.out = 101)
  first <- density_level_set_confidence(
    x, 0.25, points, bandwidth = 0.35, n_boot = 15, random_state = 8
  )
  second <- density_level_set_confidence(
    x, 0.25, points, bandwidth = 0.35, n_boot = 15, random_state = 8
  )
  expect_gte(first$radius, 0)
  expect_identical(first$bootstrap_statistics, second$bootstrap_statistics)
})

test_that("inverse regression confidence can invert a band", {
  x <- seq(-1, 1, length.out = 35)
  result <- inverse_regression_confidence(
    x, x, 0, seq(-0.8, 0.8, length.out = 51), bandwidth = 0.4,
    n_boot = 10, method = "inversion", random_state = 2
  )
  expect_true(any(result$mask))
  expect_equal(result$roots[1], 0, tolerance = 1e-8)
})

test_that("inverse regression supports Hausdorff confidence", {
  x <- seq(-1, 1, length.out = 41)
  y <- x + 0.03 * sin(8 * x)
  result <- inverse_regression_confidence(
    x, y, 0, seq(-0.8, 0.8, length.out = 61), bandwidth = 0.4,
    n_boot = 8, method = "hausdorff", random_state = 5
  )
  expect_gte(result$radius, 0)
  expect_equal(result$n_boot, 8)
  expect_length(result$bootstrap_statistics, 8)
})

test_that("inverse regression supports a normal bootstrap interval", {
  x <- seq(-1, 1, length.out = 41)
  y <- x + 0.03 * sin(8 * x)
  result <- inverse_regression_confidence(
    x, y, 0, seq(-0.8, 0.8, length.out = 61), bandwidth = 0.4,
    n_boot = 8, method = "normal", random_state = 5
  )
  expect_identical(result$method, "normal")
  expect_gte(result$radius, 0)
})

test_that("density level sets support band inversion", {
  x <- seq(-1, 1, length.out = 31)
  result <- density_level_set_confidence(
    x, 0.2, seq(-2, 2, length.out = 81), bandwidth = 0.4,
    n_boot = 8, method = "inversion", random_state = 6
  )
  expect_identical(result$method, "band_inversion")
  expect_equal(result$confidence, 0.95)
  expect_length(result$bootstrap_statistics, 8)
})
