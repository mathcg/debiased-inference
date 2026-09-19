test_that("debiased KDE matches the closed form at the origin", {
  x <- c(-1, 0, 1)
  result <- debiased_kde(x, 0, bandwidth = 1, tau = 1)
  phi <- dnorm(x)
  expected <- mean((1.5 - 0.5 * x^2) * phi)
  expect_equal(result$estimate, expected, tolerance = 1e-14)
})

test_that("debiased Gaussian kernel has unit mass and zero second moment", {
  x <- c(0, 1)
  grid <- seq(-8, 8, length.out = 20001)
  estimate <- debiased_kde(x, grid, bandwidth = 1, tau = 1)$estimate
  dx <- grid[2] - grid[1]
  expect_equal(sum(estimate) * dx, 1, tolerance = 1e-9)
  expect_equal(sum(grid^2 * estimate) * dx, mean(x^2), tolerance = 1e-8)
})

test_that("KDE matches the shared cross-language reference", {
  result <- debiased_kde(
    c(-1, 0, 1), c(-0.5, 0, 0.5), bandwidth = 0.7, tau = 0.8
  )
  expected <- c(
    0.34746560868427934,
    0.36005927354800898,
    0.34746560868427934
  )
  expect_equal(result$estimate, expected, tolerance = 1e-14)
})

test_that("multivariate KDE returns finite values", {
  x <- rbind(c(0, 0), c(1, -1), c(-0.5, 0.5))
  points <- rbind(c(0, 0), c(0.2, -0.1))
  result <- debiased_kde(x, points, bandwidth = 0.8)
  expect_equal(dim(result$points), c(2, 2))
  expect_length(result$estimate, 2)
  expect_true(all(is.finite(result$estimate)))
})

test_that("density CV bandwidth selects a candidate", {
  x <- c(-1.2, -0.9, -0.4, 0.1, 0.2, 0.8, 1.4)
  candidates <- c(0.2, 0.4, 0.8)
  selected <- density_bandwidth(x, method = "cv", candidates = candidates)
  expect_true(selected %in% candidates)
})

test_that("fixed-width bootstrap band is reproducible and symmetric", {
  x <- c(-1.2, -0.7, -0.2, 0.1, 0.4, 0.9, 1.3)
  points <- seq(-1, 1, length.out = 15)
  first <- kde_confidence_band(
    x, points, bandwidth = 0.5, n_boot = 30, random_state = 42
  )
  second <- kde_confidence_band(
    x, points, bandwidth = 0.5, n_boot = 30, random_state = 42
  )
  expect_identical(first$bootstrap_statistics, second$bootstrap_statistics)
  expect_equal(first$upper - first$estimate,
               rep(first$critical_value, length(points)))
  expect_equal(first$estimate - first$lower,
               rep(first$critical_value, length(points)))
})

test_that("invalid bandwidths are rejected", {
  for (bad in c(0, -1, Inf, NaN)) {
    expect_error(debiased_kde(c(0, 1), 0, bandwidth = bad), "bandwidth")
  }
})

test_that("studentized must be logical", {
  expect_error(kde_confidence_band(
    c(-1, 0, 1), 0, bandwidth = 0.5, n_boot = 2,
    studentized = "yes"
  ), "TRUE or FALSE")
})
