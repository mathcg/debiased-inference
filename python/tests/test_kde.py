import numpy as np
import pytest

from debiased_inference import (
    debiased_kde,
    density_bandwidth,
    kde_confidence_band,
)


def test_debiased_kde_matches_closed_form_at_origin():
    x = np.array([-1.0, 0.0, 1.0])
    result = debiased_kde(x, [0.0], bandwidth=1.0, tau=1.0)
    phi = np.exp(-0.5 * x**2) / np.sqrt(2.0 * np.pi)
    expected = np.mean((1.5 - 0.5 * x**2) * phi)
    np.testing.assert_allclose(result.estimate, [expected], rtol=1e-14)


def test_cross_language_kde_reference_values():
    result = debiased_kde(
        [-1.0, 0.0, 1.0], [-0.5, 0.0, 0.5], bandwidth=0.7, tau=0.8
    )
    expected = [
        0.34746560868427934,
        0.36005927354800898,
        0.34746560868427934,
    ]
    np.testing.assert_allclose(result.estimate, expected, rtol=1e-14)


def test_debiased_gaussian_kernel_has_unit_mass_and_zero_second_moment():
    x = np.array([0.0, 1.0])
    grid = np.linspace(-8.0, 8.0, 20001)
    result = debiased_kde(x, grid, bandwidth=1.0, tau=1.0)
    integral = np.trapezoid(result.estimate, grid)
    second_moment = np.trapezoid(grid**2 * result.estimate, grid)
    assert integral == pytest.approx(1.0, abs=1e-10)
    # M_tau has zero second moment, so smoothing does not add to E[X^2].
    assert second_moment == pytest.approx(np.mean(x**2), abs=2e-9)


def test_multivariate_estimate_shape_and_finiteness():
    x = np.array([[0.0, 0.0], [1.0, -1.0], [-0.5, 0.5]])
    points = np.array([[0.0, 0.0], [0.2, -0.1]])
    result = debiased_kde(x, points, bandwidth=0.8)
    assert result.points.shape == (2, 2)
    assert result.estimate.shape == (2,)
    assert np.all(np.isfinite(result.estimate))


def test_density_bandwidth_is_positive_for_partially_degenerate_multivariate_data():
    x = np.column_stack([np.arange(10.0), np.ones(10)])
    assert density_bandwidth(x) > 0


def test_density_cv_bandwidth_selects_a_candidate():
    x = np.array([-1.2, -0.9, -0.4, 0.1, 0.2, 0.8, 1.4])
    candidates = np.array([0.2, 0.4, 0.8])
    selected = density_bandwidth(x, method="cv", candidates=candidates)
    assert selected in candidates


def test_fixed_width_band_is_reproducible_and_symmetric():
    x = np.array([-1.2, -0.7, -0.2, 0.1, 0.4, 0.9, 1.3])
    points = np.linspace(-1.0, 1.0, 15)
    first = kde_confidence_band(
        x, points, bandwidth=0.5, n_boot=30, random_state=42
    )
    second = kde_confidence_band(
        x, points, bandwidth=0.5, n_boot=30, random_state=42
    )
    np.testing.assert_array_equal(first.bootstrap_statistics, second.bootstrap_statistics)
    np.testing.assert_allclose(first.upper - first.estimate, first.critical_value)
    np.testing.assert_allclose(first.estimate - first.lower, first.critical_value)


def test_studentized_band_has_variable_width():
    rng = np.random.default_rng(9)
    x = rng.normal(size=50)
    result = kde_confidence_band(
        x,
        np.linspace(-1.5, 1.5, 21),
        bandwidth=0.45,
        n_boot=20,
        random_state=4,
        studentized=True,
    )
    assert result.standard_error is not None
    assert np.ptp(result.width) > 0


def test_studentized_requires_boolean():
    with pytest.raises(TypeError, match="boolean"):
        kde_confidence_band(
            [-1.0, 0.0, 1.0], [0.0], bandwidth=0.5,
            n_boot=2, studentized="yes",
        )


@pytest.mark.parametrize("bad", [0, -1, np.inf, np.nan])
def test_invalid_bandwidth_is_rejected(bad):
    with pytest.raises(ValueError, match="bandwidth"):
        debiased_kde([0.0, 1.0], [0.0], bandwidth=bad)
