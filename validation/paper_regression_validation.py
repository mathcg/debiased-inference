"""Reproduce regression-band simulations from Cheng and Chen (2019).

The ``table1`` study targets the debiased-CV column in Tables 1 and 2.  The
``sine`` study targets the debiased curves in Figure 7.  Bands and bandwidths
are computed exclusively through the package's public API.
"""

from __future__ import annotations

import argparse
import csv
import math
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from debiased_inference import regression_bandwidth, regression_confidence_band

PAPER_TABLE1 = {
    500: (0.976, 0.090),
    1000: (0.976, 0.066),
    2000: (0.963, 0.049),
}


@dataclass(frozen=True)
class Task:
    study: str
    sample_size: int
    bandwidth_factors: tuple[float, ...]
    replication: int
    n_boot: int
    grid_size: int
    seed: int


def _regression_function(study: str, values: np.ndarray) -> np.ndarray:
    if study == "sine":
        return np.sin(np.pi * values)
    if study == "table1":
        denominator = 1.0 + 18.0 * values**2 * (np.sign(values) + 1.0)
        return np.sin(1.5 * np.pi * values) / denominator
    raise ValueError(f"unknown study: {study}")


def _run_one(task: Task) -> list[dict[str, float | int | str]]:
    data_seed = task.seed + 1_000_003 * task.replication
    generator = np.random.default_rng(data_seed)
    if task.study == "sine":
        x = generator.uniform(0.0, 1.0, task.sample_size)
        grid = np.linspace(0.1, 0.9, task.grid_size)
    else:
        x = generator.uniform(-1.0, 1.0, task.sample_size)
        grid = np.linspace(-0.9, 0.9, task.grid_size)
    truth = _regression_function(task.study, grid)
    y = _regression_function(task.study, x) + generator.normal(0.0, 0.1, task.sample_size)
    bandwidth = regression_bandwidth(x, y, n_folds=5, random_state=data_seed + 17)
    rows: list[dict[str, float | int | str]] = []
    for factor_index, factor in enumerate(task.bandwidth_factors):
        band = regression_confidence_band(
            x,
            y,
            grid,
            bandwidth=bandwidth * factor,
            confidence=0.95,
            n_boot=task.n_boot,
            random_state=data_seed + 7919 + factor_index,
        )
        rows.append(
            {
                "study": task.study,
                "sample_size": task.sample_size,
                "bandwidth_factor": factor,
                "replication": task.replication,
                "covered": int(np.max(np.abs(band.estimate - truth)) <= band.critical_value),
                "band_radius": band.critical_value,
                "bandwidth": bandwidth * factor,
                "sup_error": float(np.max(np.abs(band.estimate - truth))),
            }
        )
    return rows


def _wilson_interval(successes: int, total: int) -> tuple[float, float]:
    z = 1.959963984540054
    proportion = successes / total
    denominator = 1.0 + z**2 / total
    center = (proportion + z**2 / (2.0 * total)) / denominator
    half = z * math.sqrt(
        proportion * (1.0 - proportion) / total + z**2 / (4.0 * total**2)
    ) / denominator
    return center - half, center + half


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--study", choices=["table1", "sine"], default="table1")
    parser.add_argument("--sample-sizes", type=int, nargs="+", default=[500])
    parser.add_argument("--bandwidth-factors", type=float, nargs="+", default=[1.0])
    parser.add_argument("--replications", type=int, default=100)
    parser.add_argument("--n-boot", type=int, default=199)
    parser.add_argument("--grid-size", type=int, default=101)
    parser.add_argument("--seed", type=int, default=20260919)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()

    tasks = [
        Task(
            arguments.study,
            sample_size,
            tuple(arguments.bandwidth_factors),
            replication,
            arguments.n_boot,
            arguments.grid_size,
            arguments.seed + size_index * 10_000_019,
        )
        for size_index, sample_size in enumerate(arguments.sample_sizes)
        for replication in range(arguments.replications)
    ]
    if arguments.workers == 1:
        nested_rows = [_run_one(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=arguments.workers) as executor:
            nested_rows = list(executor.map(_run_one, tasks, chunksize=1))
    rows = [row for group in nested_rows for row in group]

    summaries: list[dict[str, object]] = []
    for sample_size in arguments.sample_sizes:
        for factor in arguments.bandwidth_factors:
            selected = [
                row
                for row in rows
                if row["sample_size"] == sample_size and row["bandwidth_factor"] == factor
            ]
            successes = sum(int(row["covered"]) for row in selected)
            low, high = _wilson_interval(successes, len(selected))
            summary: dict[str, object] = {
                "study": arguments.study,
                "method": "debiased",
                "sample_size": sample_size,
                "bandwidth_factor": factor,
                "replications": len(selected),
                "bootstrap_replicates": arguments.n_boot,
                "coverage": successes / len(selected),
                "coverage_ci_low": low,
                "coverage_ci_high": high,
                "mean_band_radius": float(np.mean([row["band_radius"] for row in selected])),
                "mean_bandwidth": float(np.mean([row["bandwidth"] for row in selected])),
                "paper_coverage": "",
                "paper_mean_band_radius": "",
            }
            if arguments.study == "table1" and factor == 1.0:
                summary["paper_coverage"], summary["paper_mean_band_radius"] = PAPER_TABLE1[
                    sample_size
                ]
            summaries.append(summary)
            print(
                f"{arguments.study} debiased n={sample_size} h_factor={factor}: "
                f"coverage={summary['coverage']:.3f} "
                f"95% CI=[{low:.3f}, {high:.3f}], "
                f"mean_radius={summary['mean_band_radius']:.4f}"
            )

    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        with arguments.output.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=list(summaries[0]), lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(summaries)


if __name__ == "__main__":
    main()
