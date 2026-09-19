# Reproduce Appendix D, Tables 3 and 4, Model 2 in Cheng and Chen (2019).
# Every confidence band is computed through the public package API.

library(debiasedInference)

arguments <- commandArgs(trailingOnly = TRUE)
replications <- if (length(arguments) >= 1L) as.integer(arguments[1L]) else 200L
n_boot <- if (length(arguments) >= 2L) as.integer(arguments[2L]) else 499L
grid_size <- if (length(arguments) >= 3L) as.integer(arguments[3L]) else 201L
output <- if (length(arguments) >= 4L) arguments[4L] else "model2-r.csv"
sample_sizes <- c(500L, 1000L, 2000L)
paper_coverage <- c(`500` = 0.956, `1000` = 0.949, `2000` = 0.953)
paper_radius <- c(`500` = 0.066, `1000` = 0.052, `2000` = 0.040)

normal_density <- function(points, mean, standard_deviation) {
  stats::dnorm(points, mean = mean, sd = standard_deviation)
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

grid <- seq(-2, 2, length.out = grid_size)
truth <- 0.5 * normal_density(grid, -1, 2 / 3) +
  0.5 * normal_density(grid, 1, 2 / 3)
rows <- vector("list", length(sample_sizes))

for (size_index in seq_along(sample_sizes)) {
  sample_size <- sample_sizes[size_index]
  covered <- logical(replications)
  radii <- numeric(replications)
  bandwidths <- numeric(replications)
  for (replication in seq_len(replications)) {
    data_seed <- 20260918L + size_index * 1000003L + replication
    set.seed(data_seed)
    means <- ifelse(stats::runif(sample_size) < 0.5, -1, 1)
    sample <- stats::rnorm(sample_size, means, 2 / 3)
    bandwidth <- density_bandwidth(sample)
    band <- kde_confidence_band(
      sample,
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
  rows[[size_index]] <- data.frame(
    study = "model2",
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
    paper_coverage = paper_coverage[as.character(sample_size)],
    paper_mean_band_radius = paper_radius[as.character(sample_size)]
  )
  cat(sprintf(
    "model2 debiased n=%d: coverage=%.3f 95%% CI=[%.3f, %.3f], mean_radius=%.4f\n",
    sample_size, mean(covered), interval[1L], interval[2L], mean(radii)
  ))
}

results <- do.call(rbind, rows)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(results, output, row.names = FALSE)
