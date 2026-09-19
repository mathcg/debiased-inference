import numpy as np
import pytest

from debiased_inference import (
    ConfidenceBandResult,
    EstimateResult,
    density_level_set_confidence,
    hausdorff_distance,
    inverse_regression_confidence,
    invert_confidence_band,
    level_set,
)
from debiased_inference.sets import _contour_points


def test_level_set_interpolates_crossings():
    estimate = EstimateResult(
        points=np.array([-1.0, 0.0, 1.0]),
        estimate=np.array([1.0, -1.0, 1.0]),
        bandwidth=0.2,
        tau=1.0,
        method="test",
    )
    result = level_set(estimate, 0.0)
    np.testing.assert_allclose(result.roots, [-0.5, 0.5])


def test_band_inversion_marks_grid_points_containing_level():
    points = np.array([-1.0, 0.0, 1.0])
    estimate = np.array([1.0, 0.0, 1.0])
    band = ConfidenceBandResult(
        points=points,
        estimate=estimate,
        lower=estimate - 0.2,
        upper=estimate + 0.2,
        critical_value=0.2,
        bandwidth=0.3,
        tau=1.0,
        confidence=0.95,
        n_boot=10,
        bootstrap_statistics=np.arange(10.0),
        method="test",
    )
    result = invert_confidence_band(band, 0.0)
    np.testing.assert_array_equal(result.mask, [False, True, False])


def test_hausdorff_distance():
    assert hausdorff_distance([0.0, 2.0], [0.5, 2.5]) == pytest.approx(0.5)


def test_two_dimensional_contour_and_hausdorff_distance():
    axis = np.linspace(-1.0, 1.0, 9)
    xx, yy = np.meshgrid(axis, axis, indexing="ij")
    points = np.column_stack([xx.ravel(), yy.ravel()])
    values = points[:, 0] ** 2 + points[:, 1] ** 2
    contour = _contour_points(points, values, 0.5)
    assert contour.shape[1] == 2
    assert contour.shape[0] > 8
    assert hausdorff_distance(contour, contour) == pytest.approx(0.0)


def test_hausdorff_rejects_empty_sets():
    with pytest.raises(ValueError):
        hausdorff_distance([], [1.0])


def test_density_level_set_hausdorff_confidence_is_reproducible():
    x = np.r_[np.linspace(-1.2, -0.3, 20), np.linspace(0.3, 1.2, 20)]
    points = np.linspace(-2.0, 2.0, 101)
    first = density_level_set_confidence(
        x,
        0.25,
        points,
        bandwidth=0.35,
        n_boot=15,
        random_state=8,
    )
    second = density_level_set_confidence(
        x,
        0.25,
        points,
        bandwidth=0.35,
        n_boot=15,
        random_state=8,
    )
    assert first.radius is not None and first.radius >= 0
    np.testing.assert_array_equal(first.bootstrap_statistics, second.bootstrap_statistics)


def test_two_dimensional_density_level_set_confidence():
    rng = np.random.default_rng(7)
    sample = rng.normal(size=(60, 2))
    axis = np.linspace(-2.5, 2.5, 17)
    xx, yy = np.meshgrid(axis, axis, indexing="ij")
    points = np.column_stack([xx.ravel(), yy.ravel()])
    result = density_level_set_confidence(
        sample,
        0.05,
        points,
        bandwidth=0.7,
        n_boot=5,
        random_state=4,
    )
    assert result.geometry.ndim == 2
    assert result.geometry.shape[1] == 2
    assert result.mask.shape == (points.shape[0],)


def test_inverse_regression_confidence_via_band_inversion():
    x = np.linspace(-1.0, 1.0, 35)
    y = x
    result = inverse_regression_confidence(
        x,
        y,
        0.0,
        np.linspace(-0.8, 0.8, 51),
        bandwidth=0.4,
        n_boot=10,
        method="inversion",
        random_state=2,
    )
    assert result.mask.any()
    assert result.roots[0] == pytest.approx(0.0, abs=1e-9)


def test_inverse_regression_hausdorff_confidence():
    x = np.linspace(-1.0, 1.0, 41)
    y = x + 0.03 * np.sin(8.0 * x)
    result = inverse_regression_confidence(
        x,
        y,
        0.0,
        np.linspace(-0.8, 0.8, 61),
        bandwidth=0.4,
        n_boot=8,
        method="hausdorff",
        random_state=5,
    )
    assert result.radius is not None and result.radius >= 0
    assert result.n_boot == 8
    assert result.bootstrap_statistics is not None


def test_inverse_regression_normal_confidence():
    x = np.linspace(-1.0, 1.0, 41)
    y = x + 0.03 * np.sin(8.0 * x)
    result = inverse_regression_confidence(
        x,
        y,
        0.0,
        np.linspace(-0.8, 0.8, 61),
        bandwidth=0.4,
        n_boot=8,
        method="normal",
        random_state=5,
    )
    assert result.method == "normal"
    assert result.radius is not None and result.radius >= 0


def test_density_level_set_confidence_via_band_inversion():
    x = np.linspace(-1.0, 1.0, 31)
    result = density_level_set_confidence(
        x,
        0.2,
        np.linspace(-2.0, 2.0, 81),
        bandwidth=0.4,
        n_boot=8,
        method="inversion",
        random_state=6,
    )
    assert result.method == "band_inversion"
    assert result.confidence == pytest.approx(0.95)
    assert result.bootstrap_statistics is not None
