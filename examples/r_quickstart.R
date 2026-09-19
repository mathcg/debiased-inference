# Reproducible R example for both confidence-band procedures.

library(debiasedInference)

set.seed(2026)
density_sample <- c(rnorm(180, -1, 0.55), rnorm(120, 1, 0.4))
density_grid <- seq(-3, 3, length.out = 250)
density_band <- kde_confidence_band(
  density_sample,
  density_grid,
  n_boot = 499,
  confidence = 0.95,
  random_state = 2026
)
print(density_band)

x <- sort(runif(300, -1, 1))
y <- sin(pi * x) + rnorm(length(x), 0, 0.25)
regression_grid <- seq(-0.9, 0.9, length.out = 200)
regression_band <- regression_confidence_band(
  x,
  y,
  regression_grid,
  n_boot = 499,
  confidence = 0.95,
  random_state = 2026
)
print(regression_band)
