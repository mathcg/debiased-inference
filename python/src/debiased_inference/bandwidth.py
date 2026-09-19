"""Bandwidth selectors for ordinary KDE and local-linear regression."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from ._validation import as_samples, as_vector, positive_scalar


def density_bandwidth(
    x: ArrayLike,
    *,
    method: str = "normal_reference",
    candidates: ArrayLike | None = None,
    block_size: int = 512,
) -> float:
    """Select an isotropic bandwidth for the ordinary KDE.

    ``method="normal_reference"`` uses the Gaussian normal-scale rule used by
    the paper's ``ks`` implementation. In one dimension this is
    ``(4 / (3 n))^(1/5) * sd``. In higher dimensions the scalar scale is the
    geometric mean of positive marginal standard deviations and the
    normal-reference dimension adjustment is used.

    ``method="cv"`` minimizes the least-squares cross-validation criterion
    over ``candidates``. If candidates are omitted, a geometric grid centered
    on the normal-reference choice is used. Pairwise calculations are blocked
    to use ``O(n * block_size)`` rather than ``O(n^2)`` memory.
    """
    samples = as_samples(x)
    if method not in {"normal_reference", "cv"}:
        raise ValueError("method must be 'normal_reference' or 'cv'")
    n, dimension = samples.shape
    standard_deviation = np.std(samples, axis=0, ddof=1)
    positive = standard_deviation[standard_deviation > 0]
    if positive.size == 0:
        raise ValueError("cannot select bandwidth from zero-scale data")
    scale = float(np.exp(np.mean(np.log(positive))))
    factor = (4.0 / (dimension + 2.0)) ** (1.0 / (dimension + 4.0))
    factor *= n ** (-1.0 / (dimension + 4.0))
    reference = positive_scalar(scale * factor, name="selected bandwidth")
    if method == "normal_reference":
        if candidates is not None:
            raise ValueError("candidates are only used when method='cv'")
        return reference

    if candidates is None:
        candidate_values = reference * np.geomspace(0.35, 2.5, 31)
    else:
        candidate_values = as_vector(candidates, name="candidates")
        if np.any(candidate_values <= 0):
            raise ValueError("candidates must be positive")
    if not isinstance(block_size, (int, np.integer)) or block_size < 1:
        raise ValueError("block_size must be a positive integer")
    scores = np.empty(candidate_values.size, dtype=float)
    for index, candidate in enumerate(candidate_values):
        h = float(candidate)
        integrated_normalizer = (4.0 * np.pi * h**2) ** (-0.5 * dimension)
        ordinary_normalizer = (2.0 * np.pi * h**2) ** (-0.5 * dimension)
        integrated_sum = 0.0
        ordinary_sum = 0.0
        for start in range(0, n, int(block_size)):
            stop = min(start + int(block_size), n)
            differences = samples[start:stop, None, :] - samples[None, :, :]
            distances_squared = np.sum(differences**2, axis=-1)
            integrated_sum += integrated_normalizer * float(
                np.sum(np.exp(-distances_squared / (4.0 * h**2)))
            )
            ordinary_sum += ordinary_normalizer * float(
                np.sum(np.exp(-distances_squared / (2.0 * h**2)))
            )
        ordinary_off_diagonal = ordinary_sum - n * ordinary_normalizer
        scores[index] = (
            integrated_sum / n**2
            - 2.0 * ordinary_off_diagonal / (n * (n - 1))
        )
    return float(candidate_values[int(np.argmin(scores))])


def regression_bandwidth(
    x: ArrayLike,
    y: ArrayLike,
    *,
    candidates: ArrayLike | None = None,
    n_folds: int = 5,
    random_state: int | None = 0,
) -> float:
    """Select a local-linear bandwidth by deterministic K-fold CV.

    Candidate loss is mean squared prediction error from the ordinary
    local-linear smoother. Failed boundary fits are excluded; a candidate is
    eligible only if every held-out point can be predicted.
    """
    from .regression import _local_polynomial

    x_values = as_vector(x, name="x")
    y_values = as_vector(y, name="y")
    if x_values.size != y_values.size:
        raise ValueError("x and y must have the same length")
    if np.ptp(x_values) == 0:
        raise ValueError("x must contain at least two distinct values")
    n = x_values.size
    if not isinstance(n_folds, (int, np.integer)) or not 2 <= n_folds <= n:
        raise ValueError("n_folds must be an integer between 2 and len(x)")

    if candidates is None:
        scale = min(float(np.std(x_values, ddof=1)), float(np.ptp(x_values)) / 4.0)
        base = max(scale * n ** (-1.0 / 5.0), np.finfo(float).eps)
        candidate_values = base * np.geomspace(0.35, 2.5, 21)
    else:
        candidate_values = as_vector(candidates, name="candidates")
        if np.any(candidate_values <= 0):
            raise ValueError("candidates must be positive")

    generator = np.random.default_rng(random_state)
    order = generator.permutation(n)
    folds = np.array_split(order, n_folds)
    losses = np.full(candidate_values.size, np.inf)
    for candidate_index, candidate in enumerate(candidate_values):
        squared_errors: list[np.ndarray] = []
        valid = True
        for held_out in folds:
            keep = np.ones(n, dtype=bool)
            keep[held_out] = False
            predicted = _local_polynomial(
                x_values[keep], y_values[keep], x_values[held_out],
                float(candidate), degree=1, derivative=0,
            )
            if np.any(~np.isfinite(predicted)):
                valid = False
                break
            squared_errors.append((y_values[held_out] - predicted) ** 2)
        if valid:
            losses[candidate_index] = float(np.mean(np.concatenate(squared_errors)))
    if not np.any(np.isfinite(losses)):
        raise RuntimeError("all candidate bandwidths produced singular local fits")
    return float(candidate_values[int(np.argmin(losses))])
