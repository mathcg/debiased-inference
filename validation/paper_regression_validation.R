# Reproduce the debiased-CV regression result in Tables 1 and 2 of
# Cheng and Chen (2019), using only the public R package API.

library(debiasedInference)

arguments <- commandArgs(trailingOnly = TRUE)
replications <- if (length(arguments) >= 1L) as.integer(arguments[1L]) else 100L
n_boot <- if (length(arguments) >= 2L) as.integer(arguments[2L]) else 199L
grid_size <- if (length(arguments) >= 3L) as.integer(arguments[3L]) else 101L
output <- if (length(arguments) >= 4L) arguments[4L] else "table1-regression-r.csv"
sample_size <- 500L

regression_function <- function(x) {
  sin(1.5 * pi * x) / (1 + 18 * x^2 * (sign(x) + 1))
}

wilson_interval <- function(successes, total) {
  z <- stats::qnorm(0.975)
  proportion <- successes / total
  denominator <- 1 + z^2 / total
  center <- (proportion + z^2 / (2 * total)) / denominator
  half <- z * sqrt(
    proportion * (1 - proportion) / total + z^2 / (4 * total^2)
  ) / denominator
  c(center - half, center + half)
}

grid <- seq(-0.9, 0.9, length.out = grid_size)
truth <- regression_function(grid)
covered <- logical(replications)
radii <- numeric(replications)
bandwidths <- numeric(replications)

for (replication in seq_len(replications)) {
  data_seed <- 20260919L + replication
  set.seed(data_seed)
  x <- stats::runif(sample_size, -1, 1)
  y <- regression_function(x) + stats::rnorm(sample_size, 0, 0.1)
  bandwidth <- regression_bandwidth(x, y, n_folds = 5L, random_state = data_seed + 17L)
  band <- regression_confidence_band(
    x,
    y,
    grid,
    bandwidth = bandwidth,
    confidence = 0.95,
    n_boot = n_boot,
    random_state = data_seed + 7919L
  )
  covered[replication] <- max(abs(band$estimate - truth)) <= band$critical_value
  radii[replication] <- band$critical_value
  bandwidths[replication] <- bandwidth
}

interval <- wilson_interval(sum(covered), replications)
results <- data.frame(
  study = "table1",
  method = "debiased",
  sample_size = sample_size,
  bandwidth_factor = 1,
  replications = replications,
  bootstrap_replicates = n_boot,
  coverage = mean(covered),
  coverage_ci_low = interval[1L],
  coverage_ci_high = interval[2L],
  mean_band_radius = mean(radii),
  mean_bandwidth = mean(bandwidths),
  paper_coverage = 0.976,
  paper_mean_band_radius = 0.090
)
cat(sprintf(
  "table1 debiased n=%d: coverage=%.3f 95%% CI=[%.3f, %.3f], mean_radius=%.4f\n",
  sample_size, mean(covered), interval[1L], interval[2L], mean(radii)
))
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(results, output, row.names = FALSE)
