"""Shared input validation helpers."""

from __future__ import annotations

from typing import Optional

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]


def as_samples(x: ArrayLike, *, name: str = "x") -> FloatArray:
    """Return observations as an ``(n, d)`` float array."""
    values = np.asarray(x, dtype=float)
    if values.ndim == 1:
        values = values[:, None]
    if values.ndim != 2 or values.shape[0] < 2 or values.shape[1] < 1:
        raise ValueError(f"{name} must contain at least two observations")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{name} must contain only finite values")
    if np.all(np.ptp(values, axis=0) == 0):
        raise ValueError(f"{name} must contain at least two distinct observations")
    return values


def as_vector(x: ArrayLike, *, name: str) -> FloatArray:
    """Return a finite, non-empty one-dimensional float array."""
    values = np.asarray(x, dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise ValueError(f"{name} must be a non-empty one-dimensional array")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{name} must contain only finite values")
    return values


def as_evaluation_points(
    points: Optional[ArrayLike], samples: FloatArray, *, grid_size: int = 200
) -> FloatArray:
    """Validate or generate evaluation points with the samples' dimension."""
    dimension = samples.shape[1]
    if points is None:
        if dimension != 1:
            raise ValueError("evaluation points are required for multivariate data")
        if not isinstance(grid_size, (int, np.integer)) or grid_size < 2:
            raise ValueError("grid_size must be an integer of at least 2")
        spread = float(np.std(samples[:, 0], ddof=1))
        pad = 0.05 * max(float(np.ptp(samples[:, 0])), spread)
        if pad == 0:
            pad = 1.0
        return np.linspace(samples[:, 0].min() - pad, samples[:, 0].max() + pad, grid_size)[:, None]

    result = np.asarray(points, dtype=float)
    if dimension == 1 and result.ndim == 1:
        result = result[:, None]
    if result.ndim != 2 or result.shape[1] != dimension or result.shape[0] == 0:
        raise ValueError(f"evaluation points must have shape (m, {dimension})")
    if not np.all(np.isfinite(result)):
        raise ValueError("evaluation points must contain only finite values")
    return result


def positive_scalar(value: float, *, name: str) -> float:
    """Validate a finite positive scalar."""
    result = float(value)
    if not np.isfinite(result) or result <= 0:
        raise ValueError(f"{name} must be a finite positive scalar")
    return result


def bootstrap_parameters(confidence: float, n_boot: int) -> tuple[float, int]:
    """Validate confidence level and bootstrap count."""
    confidence = float(confidence)
    if not 0 < confidence < 1:
        raise ValueError("confidence must be strictly between 0 and 1")
    if not isinstance(n_boot, (int, np.integer)) or n_boot < 1:
        raise ValueError("n_boot must be a positive integer")
    return confidence, int(n_boot)


def public_points(points: FloatArray) -> FloatArray:
    """Represent one-dimensional grids as vectors in public results."""
    return points[:, 0].copy() if points.shape[1] == 1 else points.copy()
