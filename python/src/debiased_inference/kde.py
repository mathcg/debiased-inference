"""Debiased Gaussian kernel density estimation and bootstrap bands."""

from __future__ import annotations

from typing import Optional

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._results import ConfidenceBandResult, EstimateResult
from ._validation import (
    FloatArray,
    as_evaluation_points,
    as_samples,
    bootstrap_parameters,
    positive_scalar,
    public_points,
)
from .bandwidth import density_bandwidth


def _gaussian_debiased_kernel(scaled_differences: FloatArray, tau: float) -> FloatArray:
    """Evaluate M_tau for an array whose last axis contains coordinates."""
    dimension = scaled_differences.shape[-1]
    norm_squared = np.sum(scaled_differences**2, axis=-1)
    normalizer = (2.0 * np.pi) ** (-0.5 * dimension)
    ordinary = normalizer * np.exp(-0.5 * norm_squared)
    tau_norm_squared = tau * tau * norm_squared
    laplacian_at_tau = normalizer * np.exp(-0.5 * tau_norm_squared)
    laplacian_at_tau *= tau_norm_squared - dimension
    return ordinary - 0.5 * tau ** (dimension + 2) * laplacian_at_tau


def _kernel_matrix(
    samples: FloatArray, points: FloatArray, bandwidth: float, tau: float
) -> FloatArray:
    differences = (points[:, None, :] - samples[None, :, :]) / bandwidth
    return _gaussian_debiased_kernel(differences, tau)


def _estimate_from_matrix(matrix: FloatArray, bandwidth: float, dimension: int) -> FloatArray:
    return np.mean(matrix, axis=1) / bandwidth**dimension


def debiased_kde(
    x: ArrayLike,
    points: Optional[ArrayLike] = None,
    *,
    bandwidth: Optional[float] = None,
    tau: float = 1.0,
    grid_size: int = 200,
) -> EstimateResult:
    """Evaluate the paper's debiased Gaussian KDE.

    Parameters
    ----------
    x:
        One-dimensional observations or an ``(n, d)`` sample matrix.
    points:
        Evaluation vector/matrix. Required for multivariate samples.
    bandwidth:
        Positive isotropic bandwidth selected for the ordinary KDE. If omitted,
        use :func:`density_bandwidth`.
    tau:
        Ratio ``h / b``. The paper recommends the default of one.
    grid_size:
        Generated grid size for one-dimensional data when ``points`` is absent.
    """
    samples = as_samples(x)
    evaluation = as_evaluation_points(points, samples, grid_size=grid_size)
    selected = (
        density_bandwidth(samples)
        if bandwidth is None
        else positive_scalar(bandwidth, name="bandwidth")
    )
    tau = positive_scalar(tau, name="tau")
    matrix = _kernel_matrix(samples, evaluation, selected, tau)
    estimate = _estimate_from_matrix(matrix, selected, samples.shape[1])
    return EstimateResult(public_points(evaluation), estimate, selected, tau, "debiased_kde")


def kde_confidence_band(
    x: ArrayLike,
    points: Optional[ArrayLike] = None,
    *,
    bandwidth: Optional[float] = None,
    tau: float = 1.0,
    confidence: float = 0.95,
    n_boot: int = 999,
    studentized: bool = False,
    random_state: Optional[int] = None,
    grid_size: int = 200,
) -> ConfidenceBandResult:
    """Construct a simultaneous empirical-bootstrap band for a density.

    The fixed-width procedure is Figure 2 of the paper. ``studentized=True``
    implements the variable-width band in Remark 1. The original bandwidth is
    held fixed in every bootstrap resample.
    """
    confidence, n_boot = bootstrap_parameters(confidence, n_boot)
    if not isinstance(studentized, (bool, np.bool_)):
        raise ValueError("studentized must be boolean")
    samples = as_samples(x)
    evaluation = as_evaluation_points(points, samples, grid_size=grid_size)
    selected = (
        density_bandwidth(samples)
        if bandwidth is None
        else positive_scalar(bandwidth, name="bandwidth")
    )
    tau = positive_scalar(tau, name="tau")
    dimension = samples.shape[1]
    matrix = _kernel_matrix(samples, evaluation, selected, tau)
    estimate = _estimate_from_matrix(matrix, selected, dimension)

    # sigma_rbc from Remark 1 is on the sqrt(n h^d) scale. Dividing by
    # sqrt(n h^d) yields the pointwise standard error used for band widths.
    squared_matrix = matrix**2
    second_moment = np.mean(squared_matrix, axis=1)
    first_moment = np.mean(matrix, axis=1)
    sigma_squared = np.maximum((second_moment - first_moment**2) / selected**dimension, 0.0)
    standard_error = np.sqrt(sigma_squared / (samples.shape[0] * selected**dimension))
    if studentized and np.any(standard_error <= np.finfo(float).eps):
        raise ValueError("studentized band is undefined where estimated variance is zero")

    generator = np.random.default_rng(random_state)
    statistics = np.empty(n_boot, dtype=float)
    probabilities = np.full(samples.shape[0], 1.0 / samples.shape[0])
    for bootstrap_index in range(n_boot):
        counts = generator.multinomial(samples.shape[0], probabilities)
        bootstrap_estimate = (
            matrix @ counts / (samples.shape[0] * selected**dimension)
        )
        difference = np.abs(bootstrap_estimate - estimate)
        if studentized:
            boot_second = squared_matrix @ counts / samples.shape[0]
            boot_first = matrix @ counts / samples.shape[0]
            boot_sigma_squared = np.maximum(
                (boot_second - boot_first**2) / selected**dimension, 0.0
            )
            bootstrap_se = np.sqrt(
                boot_sigma_squared / (samples.shape[0] * selected**dimension)
            )
            if np.any(bootstrap_se <= np.finfo(float).eps):
                raise RuntimeError("a bootstrap resample has zero estimated variance")
            difference = difference / bootstrap_se
        statistics[bootstrap_index] = float(np.max(difference))

    critical = float(np.quantile(statistics, confidence))
    half_width: NDArray[np.float64]
    half_width = critical * standard_error if studentized else np.full(estimate.shape, critical)
    return ConfidenceBandResult(
        points=public_points(evaluation),
        estimate=estimate,
        lower=estimate - half_width,
        upper=estimate + half_width,
        critical_value=critical,
        bandwidth=selected,
        tau=tau,
        confidence=confidence,
        n_boot=n_boot,
        bootstrap_statistics=statistics,
        method="debiased_kde",
        studentized=studentized,
        standard_error=standard_error if studentized else None,
    )
