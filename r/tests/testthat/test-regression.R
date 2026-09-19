test_that("debiased local-linear smoother reproduces a line", {
  x <- seq(-1, 1, length.out = 25)
  y <- 2.5 - 1.2 * x
  points <- seq(-0.9, 0.9, length.out = 9)
  result <- debiased_local_linear(x, y, points, bandwidth = 0.35)
  expect_equal(result$estimate, 2.5 - 1.2 * points, tolerance = 1e-9)
})

test_that("regression matches the shared cross-language reference", {
  x <- seq(-1, 1, length.out = 11)
  y <- 1 + x + x^2
  result <- debiased_local_linear(
    x, y, c(-0.4, 0, 0.6), bandwidth = 0.5, tau = 1
  )
  expected <- c(
    0.67351971126373411,
    0.96040877700627969,
    1.81202982185946437
  )
  expect_equal(result$estimate, expected, tolerance = 1e-13)
})

test_that("regression bandwidth is reproducible and is a candidate", {
  x <- seq(-1, 1, length.out = 40)
  y <- sin(pi * x) + 0.02 * cos(9 * x)
  candidates <- c(0.15, 0.25, 0.4)
  first <- regression_bandwidth(
    x, y, candidates = candidates, n_folds = 4, random_state = 11
  )
  second <- regression_bandwidth(
    x, y, candidates = candidates, n_folds = 4, random_state = 11
  )
  expect_true(first %in% candidates)
  expect_identical(first, second)
})

test_that("default bandwidth search is not truncated for a curved signal", {
  set.seed(7)
  x <- stats::runif(500, -1, 1)
  denominator <- 1 + 18 * x^2 * (sign(x) + 1)
  y <- sin(1.5 * pi * x) / denominator + stats::rnorm(500, 0, 0.1)
  base <- min(stats::sd(x), diff(range(x)) / 4) * length(x)^(-0.2)
  selected <- regression_bandwidth(x, y, random_state = 8)
  expect_lt(selected, 0.35 * base)
})

test_that("paired-bootstrap regression band is fixed-width", {
  x <- seq(-1, 1, length.out = 35)
  y <- sin(pi * x) + 0.1 * cos(7 * x)
  points <- seq(-0.8, 0.8, length.out = 11)
  result <- regression_confidence_band(
    x, y, points, bandwidth = 0.4, n_boot = 20, random_state = 19
  )
  expect_equal(result$upper - result$lower,
               rep(2 * result$critical_value, length(points)))
})
