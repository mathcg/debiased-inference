"""Small, configurable reproduction of the paper's density-band experiment.

The paper uses 1,000 Monte Carlo replications and 1,000 bootstrap samples.
Defaults here are intentionally smaller for an interactive smoke run; pass
``--replications 1000 --n-boot 1000`` for the paper-scale experiment.
"""

from __future__ import annotations

import argparse

import numpy as np

from debiased_inference import density_bandwidth, kde_confidence_band


def normal_density(x: np.ndarray, mean: float) -> np.ndarray:
    return np.exp(-0.5 * (x - mean) ** 2) / np.sqrt(2.0 * np.pi)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-size", type=int, default=500)
    parser.add_argument("--replications", type=int, default=10)
    parser.add_argument("--n-boot", type=int, default=199)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--bandwidth-factor", type=float, default=1.0)
    arguments = parser.parse_args()

    generator = np.random.default_rng(arguments.seed)
    grid = np.linspace(-4.0, 8.0, 301)
    truth = 0.6 * normal_density(grid, 0.0) + 0.4 * normal_density(grid, 4.0)
    covered = 0
    widths: list[float] = []
    for replication in range(arguments.replications):
        component = generator.random(arguments.sample_size) >= 0.6
        sample = generator.normal(size=arguments.sample_size) + 4.0 * component
        bandwidth = density_bandwidth(sample) * arguments.bandwidth_factor
        band = kde_confidence_band(
            sample,
            grid,
            bandwidth=bandwidth,
            n_boot=arguments.n_boot,
            random_state=arguments.seed + replication + 1,
        )
        covered += int(np.all((band.lower <= truth) & (truth <= band.upper)))
        widths.append(float(np.mean(band.width)))

    print(f"coverage={covered / arguments.replications:.4f}")
    print(f"mean_band_width={np.mean(widths):.6f}")


if __name__ == "__main__":
    main()
