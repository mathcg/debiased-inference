"""Grid-based density level sets and inverse-regression sets."""

from __future__ import annotations

from statistics import NormalDist

import numpy as np
from numpy.typing import ArrayLike

from ._results import ConfidenceBandResult, EstimateResult, SetEstimateResult
from ._validation import (
    as_evaluation_points,
    as_samples,
    as_vector,
    bootstrap_parameters,
    positive_scalar,
)


def _crossing_roots(points: np.ndarray, values: np.ndarray, level: float) -> np.ndarray:
    order = np.argsort(points)
    x = points[order]
    shifted = values[order] - level
    roots: list[float] = []
    exact = np.flatnonzero(shifted == 0)
    roots.extend(float(x[index]) for index in exact)
    crossing = np.flatnonzero(shifted[:-1] * shifted[1:] < 0)
    for index in crossing:
        fraction = -shifted[index] / (shifted[index + 1] - shifted[index])
        roots.append(float(x[index] + fraction * (x[index + 1] - x[index])))
    return np.unique(np.asarray(roots, dtype=float))


def _contour_points(
    points: np.ndarray, values: np.ndarray, level: float
) -> np.ndarray:
    """Approximate a 2-D contour on a complete rectangular grid."""
    x_coordinates = np.unique(points[:, 0])
    y_coordinates = np.unique(points[:, 1])
    if points.shape[0] != x_coordinates.size * y_coordinates.size:
        raise ValueError("two-dimensional level sets require a complete rectangular grid")
    surface = np.full((x_coordinates.size, y_coordinates.size), np.nan)
    x_index = np.searchsorted(x_coordinates, points[:, 0])
    y_index = np.searchsorted(y_coordinates, points[:, 1])
    surface[x_index, y_index] = values
    if np.any(~np.isfinite(surface)):
        raise ValueError("two-dimensional level sets require unique rectangular grid points")

    intersections: list[np.ndarray] = []
    for i in range(x_coordinates.size - 1):
        for j in range(y_coordinates.size - 1):
            vertices = np.array(
                [
                    [x_coordinates[i], y_coordinates[j]],
                    [x_coordinates[i + 1], y_coordinates[j]],
                    [x_coordinates[i + 1], y_coordinates[j + 1]],
                    [x_coordinates[i], y_coordinates[j + 1]],
                ]
            )
            heights = np.array(
                [
                    surface[i, j],
                    surface[i + 1, j],
                    surface[i + 1, j + 1],
                    surface[i, j + 1],
                ]
            )
            for triangle in ((0, 1, 2), (0, 2, 3)):
                triangle_points = vertices[list(triangle)]
                shifted = heights[list(triangle)] - level
                found: list[np.ndarray] = []
                for left, right in ((0, 1), (1, 2), (2, 0)):
                    if shifted[left] == 0:
                        found.append(triangle_points[left])
                    if shifted[left] * shifted[right] < 0:
                        fraction = -shifted[left] / (shifted[right] - shifted[left])
                        found.append(
                            triangle_points[left]
                            + fraction
                            * (triangle_points[right] - triangle_points[left])
                        )
                if shifted[2] == 0:
                    found.append(triangle_points[2])
                if len(found) >= 2:
                    intersections.extend(found[:2])
    if not intersections:
        return np.empty((0, 2), dtype=float)
    return np.unique(np.round(np.vstack(intersections), decimals=14), axis=0)


def _level_geometry(points: np.ndarray, values: np.ndarray, level: float) -> np.ndarray:
    if points.ndim == 1:
        return _crossing_roots(points, values, level)
    if points.ndim == 2 and points.shape[1] == 2:
        return _contour_points(points, values, level)
    raise ValueError("level-set geometry is supported for one or two dimensions")


def _point_cloud(values: ArrayLike, *, name: str) -> np.ndarray:
    result = np.asarray(values, dtype=float)
    if result.ndim == 1:
        result = result[:, None]
    if result.ndim != 2 or result.shape[0] == 0 or result.shape[1] == 0:
        raise ValueError(f"{name} must be a non-empty vector or point matrix")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    return result


def _distance_to_geometry(points: np.ndarray, geometry: np.ndarray) -> np.ndarray:
    point_cloud = points[:, None] if points.ndim == 1 else points
    geometry_cloud = geometry[:, None] if geometry.ndim == 1 else geometry
    differences = point_cloud[:, None, :] - geometry_cloud[None, :, :]
    return np.min(np.sqrt(np.sum(differences**2, axis=-1)), axis=1)


def level_set(estimate: EstimateResult, level: float) -> SetEstimateResult:
    """Estimate a 1-D or regular-grid 2-D equality level set."""
    points = np.asarray(estimate.points, dtype=float)
    if (
        points.ndim not in {1, 2}
        or points.shape[0] == 0
        or not np.all(np.isfinite(points))
    ):
        raise ValueError("estimate.points must be a finite vector or point matrix")
    level = float(level)
    if not np.isfinite(level):
        raise ValueError("level must be finite")
    roots = _level_geometry(points, estimate.estimate, level)
    return SetEstimateResult(
        points=points.copy(),
        mask=np.zeros(points.size, dtype=bool),
        roots=roots,
        level=level,
        method="level_set",
    )


def invert_confidence_band(
    band: ConfidenceBandResult, level: float
) -> SetEstimateResult:
    """Invert a simultaneous band to obtain a confidence set for a level set.

    This implements the alternative construction in Remark 2 for density and
    its regression analogue in Section 3.2.1.
    """
    points = np.asarray(band.points, dtype=float)
    if (
        points.ndim not in {1, 2}
        or points.shape[0] == 0
        or not np.all(np.isfinite(points))
    ):
        raise ValueError("band.points must be a finite vector or point matrix")
    level = float(level)
    if not np.isfinite(level):
        raise ValueError("level must be finite")
    mask = (band.lower <= level) & (level <= band.upper)
    roots = _level_geometry(points, band.estimate, level)
    return SetEstimateResult(
        points=points.copy(),
        mask=mask,
        roots=roots,
        level=level,
        method="band_inversion",
        confidence=band.confidence,
        n_boot=band.n_boot,
        bootstrap_statistics=band.bootstrap_statistics.copy(),
    )


def hausdorff_distance(a: ArrayLike, b: ArrayLike) -> float:
    """Hausdorff distance between finite point clouds of equal dimension."""
    left = _point_cloud(a, name="a")
    right = _point_cloud(b, name="b")
    if left.shape[1] != right.shape[1]:
        raise ValueError("a and b must have the same point dimension")
    differences = left[:, None, :] - right[None, :, :]
    pairwise = np.sqrt(np.sum(differences**2, axis=-1))
    return float(max(np.max(np.min(pairwise, axis=1)), np.max(np.min(pairwise, axis=0))))


def density_level_set(
    x: ArrayLike,
    level: float,
    points: ArrayLike | None = None,
    *,
    bandwidth: float | None = None,
    tau: float = 1.0,
    grid_size: int = 200,
) -> SetEstimateResult:
    """Estimate a 1-D or regular-grid 2-D density equality level set."""
    from .kde import debiased_kde

    return level_set(
        debiased_kde(
            x, points, bandwidth=bandwidth, tau=tau, grid_size=grid_size
        ),
        level,
    )


def density_level_set_confidence(
    x: ArrayLike,
    level: float,
    points: ArrayLike | None = None,
    *,
    bandwidth: float | None = None,
    tau: float = 1.0,
    confidence: float = 0.95,
    n_boot: int = 999,
    method: str = "hausdorff",
    random_state: int | None = None,
    grid_size: int = 200,
    max_attempts: int | None = None,
) -> SetEstimateResult:
    """Construct a confidence set for a 1-D or regular-grid 2-D density level set.

    ``method="hausdorff"`` implements Section 3.1.1 by bootstrapping the
    Hausdorff distance between interpolated level-set roots. ``"inversion"``
    implements Remark 2 by inverting the simultaneous density band.
    """
    from .bandwidth import density_bandwidth
    from .kde import _estimate_from_matrix, _kernel_matrix, kde_confidence_band

    if method not in {"hausdorff", "inversion"}:
        raise ValueError("method must be 'hausdorff' or 'inversion'")
    if method == "inversion":
        return invert_confidence_band(
            kde_confidence_band(
                x,
                points,
                bandwidth=bandwidth,
                tau=tau,
                confidence=confidence,
                n_boot=n_boot,
                random_state=random_state,
                grid_size=grid_size,
            ),
            level,
        )

    confidence, n_boot = bootstrap_parameters(confidence, n_boot)
    samples = as_samples(x)
    if samples.shape[1] not in {1, 2}:
        raise ValueError("Hausdorff level-set confidence supports one or two dimensions")
    evaluation_matrix = as_evaluation_points(points, samples, grid_size=grid_size)
    evaluation = (
        evaluation_matrix[:, 0]
        if samples.shape[1] == 1
        else evaluation_matrix.copy()
    )
    selected = density_bandwidth(samples) if bandwidth is None else positive_scalar(
        bandwidth, name="bandwidth"
    )
    tau = positive_scalar(tau, name="tau")
    level = float(level)
    if not np.isfinite(level):
        raise ValueError("level must be finite")
    matrix = _kernel_matrix(samples, evaluation_matrix, selected, tau)
    estimate = _estimate_from_matrix(matrix, selected, samples.shape[1])
    roots = _level_geometry(evaluation, estimate, level)
    if roots.size == 0:
        raise ValueError("estimated level set is empty on the evaluation grid")

    generator = np.random.default_rng(random_state)
    attempts_limit = max_attempts if max_attempts is not None else max(10 * n_boot, 100)
    if not isinstance(attempts_limit, (int, np.integer)) or attempts_limit < n_boot:
        raise ValueError("max_attempts must be an integer at least n_boot")
    statistics = np.empty(n_boot, dtype=float)
    accepted = 0
    attempts = 0
    probabilities = np.full(samples.shape[0], 1.0 / samples.shape[0])
    while accepted < n_boot and attempts < attempts_limit:
        attempts += 1
        counts = generator.multinomial(samples.shape[0], probabilities)
        bootstrap_estimate = (
            matrix @ counts
            / (samples.shape[0] * selected ** samples.shape[1])
        )
        bootstrap_roots = _level_geometry(evaluation, bootstrap_estimate, level)
        if bootstrap_roots.size:
            statistics[accepted] = hausdorff_distance(bootstrap_roots, roots)
            accepted += 1
    if accepted < n_boot:
        raise RuntimeError(
            f"only {accepted} non-empty bootstrap level sets after {attempts} attempts"
        )
    radius = float(np.quantile(statistics, confidence))
    mask = _distance_to_geometry(evaluation, roots) <= radius
    return SetEstimateResult(
        points=evaluation.copy(),
        mask=mask,
        roots=roots,
        level=level,
        method="hausdorff",
        radius=radius,
        confidence=confidence,
        n_boot=n_boot,
        bootstrap_statistics=statistics,
    )


def inverse_regression(
    x: ArrayLike,
    y: ArrayLike,
    level: float,
    points: ArrayLike | None = None,
    *,
    bandwidth: float | None = None,
    tau: float = 1.0,
    grid_size: int = 200,
    n_folds: int = 5,
    random_state: int | None = 0,
) -> SetEstimateResult:
    """Estimate a one-dimensional inverse-regression equality set."""
    from .regression import debiased_local_linear

    return level_set(
        debiased_local_linear(
            x,
            y,
            points,
            bandwidth=bandwidth,
            tau=tau,
            grid_size=grid_size,
            n_folds=n_folds,
            random_state=random_state,
        ),
        level,
    )


def inverse_regression_confidence(
    x: ArrayLike,
    y: ArrayLike,
    level: float,
    points: ArrayLike | None = None,
    *,
    bandwidth: float | None = None,
    tau: float = 1.0,
    confidence: float = 0.95,
    n_boot: int = 999,
    method: str = "hausdorff",
    random_state: int | None = None,
    grid_size: int = 200,
    n_folds: int = 5,
    max_attempts: int | None = None,
) -> SetEstimateResult:
    """Construct a confidence set for an inverse-regression equality set."""
    from .regression import (
        _debiased_regression_values,
        debiased_local_linear,
        regression_confidence_band,
    )

    if method not in {"hausdorff", "inversion", "normal"}:
        raise ValueError("method must be 'hausdorff', 'inversion', or 'normal'")
    if method == "inversion":
        return invert_confidence_band(
            regression_confidence_band(
                x,
                y,
                points,
                bandwidth=bandwidth,
                tau=tau,
                confidence=confidence,
                n_boot=n_boot,
                random_state=random_state,
                grid_size=grid_size,
                n_folds=n_folds,
                max_attempts=max_attempts,
            ),
            level,
        )

    confidence, n_boot = bootstrap_parameters(confidence, n_boot)
    x_values = as_vector(x, name="x")
    y_values = as_vector(y, name="y")
    if x_values.size != y_values.size or x_values.size < 4:
        raise ValueError("x and y must have equal lengths of at least 4")
    base = debiased_local_linear(
        x_values,
        y_values,
        points,
        bandwidth=bandwidth,
        tau=tau,
        grid_size=grid_size,
        n_folds=n_folds,
        random_state=random_state,
    )
    level = float(level)
    if not np.isfinite(level):
        raise ValueError("level must be finite")
    roots = _crossing_roots(base.points, base.estimate, level)
    if roots.size == 0:
        raise ValueError("estimated inverse-regression set is empty on the evaluation grid")
    if method == "normal" and roots.size != 1:
        raise ValueError("method='normal' requires exactly one estimated crossing")
    if method == "normal" and n_boot < 2:
        raise ValueError("method='normal' requires at least two bootstrap replicates")

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
            bootstrap_roots = _crossing_roots(base.points, estimate, level)
            valid_roots = (
                bootstrap_roots.size == 1
                if method == "normal"
                else bootstrap_roots.size > 0
            )
            if valid_roots:
                statistics[accepted] = (
                    bootstrap_roots[0]
                    if method == "normal"
                    else hausdorff_distance(bootstrap_roots, roots)
                )
                accepted += 1
    if accepted < n_boot:
        raise RuntimeError(
            f"only {accepted} non-empty bootstrap inverse sets after {attempts} attempts"
        )
    if method == "normal":
        standard_error = float(np.std(statistics, ddof=1))
        radius = NormalDist().inv_cdf(0.5 + confidence / 2.0) * standard_error
    else:
        radius = float(np.quantile(statistics, confidence))
    mask = np.min(np.abs(base.points[:, None] - roots[None, :]), axis=1) <= radius
    return SetEstimateResult(
        points=base.points.copy(),
        mask=mask,
        roots=roots,
        level=level,
        method=method,
        radius=radius,
        confidence=confidence,
        n_boot=n_boot,
        bootstrap_statistics=statistics,
    )
