"""Reproducible Python example for both confidence-band procedures."""

import numpy as np

from debiased_inference import kde_confidence_band, regression_confidence_band


rng = np.random.default_rng(2026)

density_sample = np.r_[rng.normal(-1.0, 0.55, 180), rng.normal(1.0, 0.4, 120)]
density_grid = np.linspace(-3.0, 3.0, 250)
density_band = kde_confidence_band(
    density_sample,
    density_grid,
    n_boot=499,
    confidence=0.95,
    random_state=2026,
)
print("Density bandwidth:", density_band.bandwidth)
print("Density critical value:", density_band.critical_value)

x = np.sort(rng.uniform(-1.0, 1.0, 300))
y = np.sin(np.pi * x) + rng.normal(0.0, 0.25, x.size)
regression_grid = np.linspace(-0.9, 0.9, 200)
regression_band = regression_confidence_band(
    x,
    y,
    regression_grid,
    n_boot=499,
    confidence=0.95,
    random_state=2026,
)
print("Regression bandwidth:", regression_band.bandwidth)
print("Regression critical value:", regression_band.critical_value)
