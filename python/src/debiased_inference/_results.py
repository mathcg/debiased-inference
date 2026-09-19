"""Immutable result containers for the public API."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class EstimateResult:
    """A nonparametric estimate evaluated on a finite grid."""

    points: FloatArray
    estimate: FloatArray
    bandwidth: float
    tau: float
    method: str

    def __repr__(self) -> str:
        return (
            f"EstimateResult(method={self.method!r}, points={len(self.points)}, "
            f"bandwidth={self.bandwidth:.6g}, tau={self.tau:.6g})"
        )


@dataclass(frozen=True)
class ConfidenceBandResult:
    """A simultaneous bootstrap confidence band on a finite grid."""

    points: FloatArray
    estimate: FloatArray
    lower: FloatArray
    upper: FloatArray
    critical_value: float
    bandwidth: float
    tau: float
    confidence: float
    n_boot: int
    bootstrap_statistics: FloatArray
    method: str
    studentized: bool = False
    standard_error: FloatArray | None = None

    @property
    def width(self) -> FloatArray:
        """Pointwise width of the simultaneous band."""
        return self.upper - self.lower

    def __repr__(self) -> str:
        return (
            f"ConfidenceBandResult(method={self.method!r}, points={len(self.points)}, "
            f"confidence={self.confidence:.3g}, n_boot={self.n_boot}, "
            f"bandwidth={self.bandwidth:.6g})"
        )


@dataclass(frozen=True)
class SetEstimateResult:
    """A grid-based estimate or confidence region for a level set.

    ``roots`` is a vector for a one-dimensional level set and an ``(k, 2)``
    contour point cloud for a two-dimensional level set.
    """

    points: FloatArray
    mask: NDArray[np.bool_]
    roots: FloatArray
    level: float
    method: str
    radius: float | None = None
    confidence: float | None = None
    n_boot: int | None = None
    bootstrap_statistics: FloatArray | None = None

    @property
    def geometry(self) -> FloatArray:
        """Return the interpolated level-set geometry."""
        return self.roots

    def __repr__(self) -> str:
        size = len(self.roots)
        detail = f", radius={self.radius:.6g}" if self.radius is not None else ""
        return (
            f"SetEstimateResult(method={self.method!r}, level={self.level:.6g}, "
            f"geometry_points={size}{detail})"
        )
