# Small, configurable reproduction of the paper's density-band experiment.
# Usage: Rscript paper_density_simulation.R [sample_size] [replications]
#                                         [n_boot] [seed] [bandwidth_factor]
# Use 500 1000 1000 2026 1 for the first paper-scale sample size.

library(debiasedInference)

arguments <- commandArgs(trailingOnly = TRUE)
sample_size <- if (length(arguments) >= 1L) as.integer(arguments[1L]) else 500L
replications <- if (length(arguments) >= 2L) as.integer(arguments[2L]) else 10L
n_boot <- if (length(arguments) >= 3L) as.integer(arguments[3L]) else 199L
seed <- if (length(arguments) >= 4L) as.integer(arguments[4L]) else 2026L
bandwidth_factor <- if (length(arguments) >= 5L) as.double(arguments[5L]) else 1

set.seed(seed)
grid <- seq(-4, 8, length.out = 301L)
truth <- 0.6 * dnorm(grid) + 0.4 * dnorm(grid, mean = 4)
covered <- logical(replications)
widths <- numeric(replications)

for (replication in seq_len(replications)) {
  component <- runif(sample_size) >= 0.6
  sample <- rnorm(sample_size) + 4 * component
  bandwidth <- density_bandwidth(sample) * bandwidth_factor
  band <- kde_confidence_band(
    sample, grid, bandwidth = bandwidth, n_boot = n_boot,
    random_state = seed + replication
  )
  covered[replication] <- all(band$lower <= truth & truth <= band$upper)
  widths[replication] <- mean(band$upper - band$lower)
}

cat(sprintf("coverage=%.4f\n", mean(covered)))
cat(sprintf("mean_band_radius=%.6f\n", mean(widths) / 2))
