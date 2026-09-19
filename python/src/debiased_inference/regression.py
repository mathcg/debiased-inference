"""Debiased local-linear regression and paired-bootstrap bands."""

from __future__ import annotations

import math
from typing import Optional

import numpy as np
from numpy.typing import ArrayLike

from ._results import ConfidenceBandResult, EstimateResult
from ._validation import (
    FloatArray,
    as_vector,
    bootstrap_parameters,
    positive_scalar,
)


def _local_polynomial(
    x: FloatArray,
    y: FloatArray,
    points: FloatArray,
    bandwidth: float,
    *,
    degree: int,
    derivative: int,
) -> FloatArray:
    """Gaussian weighted local-polynomial derivative estimates."""
    estimates = np.full(points.size, np.nan, dtype=float)
    for point_index, point in enumerate(points):
        offsets = x - point
        scaled = offsets / bandwidth
        weights = np.exp(-0.5 * scaled**2) / math.sqrt(2.0 * math.pi)
        design = np.column_stack([offsets**power for power in range(degree + 1)])
        root_weights = np.sqrt(weights)
        weighted_design = design * root_weights[:, None]
        weighted_response = y * root_weights
        coefficients, _, rank, _ = np.linalg.lstsq(weighted_design, weighted_response, rcond=None)
        if rank == degree + 1:
            estimates[point_index] = math.factorial(derivative) * coefficients[derivative]
    return estimates


def _debiased_regression_values(
    x: FloatArray,
    y: FloatArray,
    points: FloatArray,
    bandwidth: float,
    tau: float,
) -> FloatArray:
    ordinary = _local_polynomial(x, y, points, bandwidth, degree=1, derivative=0)
    second = _local_polynomial(
        x, y, points, bandwidth / tau, degree=3, derivative=2
    )
    return ordinary - 0.5 * bandwidth**2 * second


def debiased_local_linear(
    x: ArrayLike,
    y: ArrayLike,
    points: Optional[ArrayLike] = None,
    *,
    bandwidth: Optional[float] = None,
    tau: float = 1.0,
    grid_size: int = 200,
    n_folds: int = 5,
    random_state: Optional[int] = 0,
) -> EstimateResult:
    """Evaluate the paper's debiased one-dimensional local-linear smoother."""
    x_values = as_vector(x, name="x")
    y_values = as_vector(y, name="y")
    if x_values.size != y_values.size or x_values.size < 4:
        raise ValueError("x and y must have equal lengths of at least 4")
    if np.ptp(x_values) == 0:
        raise ValueError("x must contain at least two distinct values")
    if points is None:
        if not isinstance(grid_size, (int, np.integer)) or grid_size < 2:
            raise ValueError("grid_size must be an integer of at least 2")
        evaluation = np.linspace(x_values.min(), x_values.max(), grid_size)
    else:
        evaluation = as_vector(points, name="points")
    if bandwidth is None:
        from .bandwidth import regression_bandwidth

        selected = regression_bandwidth(
            x_values, y_values, n_folds=n_folds, random_state=random_state
        )
    else:
        selected = positive_scalar(bandwidth, name="bandwidth")
    tau = positive_scalar(tau, name="tau")
    estimate = _debiased_regression_values(
        x_values, y_values, evaluation, selected, tau
    )
    if np.any(~np.isfinite(estimate)):
        raise RuntimeError(
            "singular local fit; increase bandwidth or restrict evaluation points"
        )
    return EstimateResult(evaluation.copy(), estimate, selected, tau, "debiased_local_linear")


def regression_confidence_band(
    x: ArrayLike,
    y: ArrayLike,
    points: Optional[ArrayLike] = None,
    *,
    bandwidth: Optional[float] = None,
    tau: float = 1.0,
    confidence: float = 0.95,
    n_boot: int = 999,
    random_state: Optional[int] = None,
    grid_size: int = 200,
    n_folds: int = 5,
    max_attempts: Optional[int] = None,
) -> ConfidenceBandResult:
    """Construct a simultaneous paired-bootstrap regression band."""
    confidence, n_boot = bootstrap_parameters(confidence, n_boot)
    base = debiased_local_linear(
        x,
        y,
        points,
        bandwidth=bandwidth,
        tau=tau,
        grid_size=grid_size,
        n_folds=n_folds,
        random_state=random_state,
    )
    x_values = as_vector(x, name="x")
    y_values = as_vector(y, name="y")
    generator = np.random.default_rng(random_state)
    attempts_limit = max_attempts if max_attempts is not None else max(10 * n_boot, 100)
    if not isinstance(attempts_limit, (int, np.integer)) or attempts_limit < n_boot:
        raise ValueError("max_attempts must be an integer at least n_boot")

    statistics = np.empty(n_boot, dtype=float)
    accepted = 0
    attempts = 0
    while accepted < n_boot and attempts < attempts_limit:
        attempts += 1
        indices = generator.integers(0, x_values.size, size=x_values.size)
        estimate = _debiased_regression_values(
            x_values[indices],
            y_values[indices],
            base.points,
            base.bandwidth,
            base.tau,
        )
        if np.all(np.isfinite(estimate)):
            statistics[accepted] = float(np.max(np.abs(estimate - base.estimate)))
            accepted += 1
    if accepted < n_boot:
        raise RuntimeError(
            f"only {accepted} nonsingular bootstrap fits after {attempts} attempts"
        )
    critical = float(np.quantile(statistics, confidence))
    return ConfidenceBandResult(
        points=base.points,
        estimate=base.estimate,
        lower=base.estimate - critical,
        upper=base.estimate + critical,
        critical_value=critical,
        bandwidth=base.bandwidth,
        tau=base.tau,
        confidence=confidence,
        n_boot=n_boot,
        bootstrap_statistics=statistics,
        method="debiased_local_linear",
    )
