"""Reproduce density-band simulations from Cheng and Chen (2019).

The default ``model2`` study targets Appendix D, Tables 3 and 4.  The ``main``
study targets Figure 4.  Every debiased band is computed through the public
``debiased_inference.kde_confidence_band`` API.
"""

from __future__ import annotations

import argparse
import csv
import math
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from debiased_inference import density_bandwidth, kde_confidence_band

PAPER_MODEL2 = {
    500: (0.956, 0.066),
    1000: (0.949, 0.052),
    2000: (0.953, 0.040),
}


@dataclass(frozen=True)
class Task:
    study: str
    sample_size: int
    bandwidth_factor: float
    method: str
    replication: int
    n_boot: int
    grid_size: int
    seed: int


def _normal_density(points: np.ndarray, mean: float, standard_deviation: float) -> np.ndarray:
    scaled = (points - mean) / standard_deviation
    return np.exp(-0.5 * scaled**2) / (standard_deviation * np.sqrt(2.0 * np.pi))


def _sample_and_truth(
    study: str, sample_size: int, grid_size: int, generator: np.random.Generator
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if study == "model2":
        grid = np.linspace(-2.0, 2.0, grid_size)
        means = np.where(generator.random(sample_size) < 0.5, -1.0, 1.0)
        sample = generator.normal(means, 2.0 / 3.0)
        truth = 0.5 * _normal_density(grid, -1.0, 2.0 / 3.0)
        truth += 0.5 * _normal_density(grid, 1.0, 2.0 / 3.0)
        return sample, grid, truth
    if study == "main":
        grid = np.linspace(-4.0, 8.0, grid_size)
        means = np.where(generator.random(sample_size) < 0.6, 0.0, 4.0)
        sample = generator.normal(means, 1.0)
        truth = 0.6 * _normal_density(grid, 0.0, 1.0)
        truth += 0.4 * _normal_density(grid, 4.0, 1.0)
        return sample, grid, truth
    raise ValueError(f"unknown study: {study}")


def _ordinary_band(
    sample: np.ndarray,
    grid: np.ndarray,
    bandwidth: float,
    confidence: float,
    n_boot: int,
    seed: int,
) -> tuple[np.ndarray, float]:
    scaled = (grid[:, None] - sample[None, :]) / bandwidth
    matrix = np.exp(-0.5 * scaled**2) / np.sqrt(2.0 * np.pi)
    estimate = np.mean(matrix, axis=1) / bandwidth
    generator = np.random.default_rng(seed)
    probabilities = np.full(sample.size, 1.0 / sample.size)
    statistics = np.empty(n_boot)
    for index in range(n_boot):
        counts = generator.multinomial(sample.size, probabilities)
        bootstrap = matrix @ counts / (sample.size * bandwidth)
        statistics[index] = np.max(np.abs(bootstrap - estimate))
    return estimate, float(np.quantile(statistics, confidence))


def _run_one(task: Task) -> dict[str, float | int | str]:
    data_seed = task.seed + 1_000_003 * task.replication
    bootstrap_seed = data_seed + 7919
    generator = np.random.default_rng(data_seed)
    sample, grid, truth = _sample_and_truth(
        task.study, task.sample_size, task.grid_size, generator
    )
    bandwidth = density_bandwidth(sample) * task.bandwidth_factor
    if task.method == "debiased":
        band = kde_confidence_band(
            sample,
            grid,
            bandwidth=bandwidth,
            confidence=0.95,
            n_boot=task.n_boot,
            random_state=bootstrap_seed,
        )
        estimate = band.estimate
        radius = band.critical_value
    else:
        estimate, radius = _ordinary_band(
            sample, grid, bandwidth, 0.95, task.n_boot, bootstrap_seed
        )
    return {
        "study": task.study,
        "method": task.method,
        "sample_size": task.sample_size,
        "bandwidth_factor": task.bandwidth_factor,
        "replication": task.replication,
        "covered": int(np.max(np.abs(estimate - truth)) <= radius),
        "band_radius": radius,
        "bandwidth": bandwidth,
        "sup_error": float(np.max(np.abs(estimate - truth))),
    }


def _wilson_interval(successes: int, total: int) -> tuple[float, float]:
    z = 1.959963984540054
    proportion = successes / total
    denominator = 1.0 + z**2 / total
    center = (proportion + z**2 / (2.0 * total)) / denominator
    half = z * math.sqrt(
        proportion * (1.0 - proportion) / total + z**2 / (4.0 * total**2)
    ) / denominator
    return center - half, center + half


def _summarize(rows: list[dict[str, float | int | str]]) -> list[dict[str, object]]:
    keys = sorted(
        {
            (
                str(row["study"]),
                str(row["method"]),
                int(row["sample_size"]),
                float(row["bandwidth_factor"]),
            )
            for row in rows
        }
    )
    summaries: list[dict[str, object]] = []
    for study, method, sample_size, factor in keys:
        selected = [
            row
            for row in rows
            if row["study"] == study
            and row["method"] == method
            and row["sample_size"] == sample_size
            and row["bandwidth_factor"] == factor
        ]
        covered = sum(int(row["covered"]) for row in selected)
        low, high = _wilson_interval(covered, len(selected))
        summary: dict[str, object] = {
            "study": study,
            "method": method,
            "sample_size": sample_size,
            "bandwidth_factor": factor,
            "replications": len(selected),
            "bootstrap_replicates": "",
            "coverage": covered / len(selected),
            "coverage_ci_low": low,
            "coverage_ci_high": high,
            "mean_band_radius": float(np.mean([row["band_radius"] for row in selected])),
            "mean_bandwidth": float(np.mean([row["bandwidth"] for row in selected])),
            "paper_coverage": "",
            "paper_mean_band_radius": "",
        }
        if study == "model2" and method == "debiased" and factor == 1.0:
            summary["paper_coverage"], summary["paper_mean_band_radius"] = PAPER_MODEL2[
                sample_size
            ]
        summaries.append(summary)
    return summaries


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--study", choices=["model2", "main"], default="model2")
    parser.add_argument("--sample-sizes", type=int, nargs="+", default=[500, 1000, 2000])
    parser.add_argument("--bandwidth-factors", type=float, nargs="+", default=[1.0])
    parser.add_argument("--methods", choices=["debiased", "ordinary"], nargs="+", default=["debiased"])
    parser.add_argument("--replications", type=int, default=200)
    parser.add_argument("--n-boot", type=int, default=499)
    parser.add_argument("--grid-size", type=int, default=201)
    parser.add_argument("--seed", type=int, default=20260918)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()

    if arguments.replications < 1 or arguments.n_boot < 1 or arguments.grid_size < 2:
        parser.error("replications and n-boot must be positive; grid-size must be at least 2")
    tasks = [
        Task(study, sample_size, factor, method, replication, arguments.n_boot,
             arguments.grid_size, arguments.seed + configuration * 10_000_019)
        for configuration, (study, sample_size, factor, method) in enumerate(
            (arguments.study, sample_size, factor, method)
            for sample_size in arguments.sample_sizes
            for factor in arguments.bandwidth_factors
            for method in arguments.methods
        )
        for replication in range(arguments.replications)
    ]
    if arguments.workers == 1:
        rows = [_run_one(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=arguments.workers) as executor:
            rows = list(executor.map(_run_one, tasks, chunksize=1))

    summaries = _summarize(rows)
    for summary in summaries:
        summary["bootstrap_replicates"] = arguments.n_boot
        print(
            f"{summary['study']} {summary['method']} n={summary['sample_size']} "
            f"h_factor={summary['bandwidth_factor']}: coverage={summary['coverage']:.3f} "
            f"95% CI=[{summary['coverage_ci_low']:.3f}, {summary['coverage_ci_high']:.3f}], "
            f"mean_radius={summary['mean_band_radius']:.4f}"
        )
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        with arguments.output.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
            writer.writeheader()
            writer.writerows(summaries)


if __name__ == "__main__":
    main()
