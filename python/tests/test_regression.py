import numpy as np
import pytest

from debiased_inference import (
    debiased_local_linear,
    regression_bandwidth,
    regression_confidence_band,
)
from debiased_inference.regression import _local_polynomial


def test_local_cubic_second_derivative_is_exact_for_cubic_polynomial():
    x = np.linspace(-2.0, 2.0, 31)
    y = 1.0 + 2.0 * x + 3.0 * x**2 - 0.5 * x**3
    points = np.array([-1.0, 0.0, 0.7])
    result = _local_polynomial(x, y, points, 0.7, degree=3, derivative=2)
    np.testing.assert_allclose(result, 6.0 - 3.0 * points, atol=1e-10)


def test_debiased_local_linear_reproduces_a_line():
    x = np.linspace(-1.0, 1.0, 25)
    y = 2.5 - 1.2 * x
    points = np.linspace(-0.9, 0.9, 9)
    result = debiased_local_linear(x, y, points, bandwidth=0.35)
    np.testing.assert_allclose(result.estimate, 2.5 - 1.2 * points, atol=1e-10)


def test_cross_language_regression_reference_values():
    x = np.linspace(-1.0, 1.0, 11)
    y = 1.0 + x + x**2
    result = debiased_local_linear(
        x, y, [-0.4, 0.0, 0.6], bandwidth=0.5, tau=1.0
    )
    expected = [
        0.67351971126373411,
        0.96040877700627969,
        1.81202982185946437,
    ]
    np.testing.assert_allclose(result.estimate, expected, rtol=1e-13)


def test_regression_bandwidth_is_a_candidate_and_reproducible():
    rng = np.random.default_rng(3)
    x = np.linspace(-1.0, 1.0, 40)
    y = np.sin(np.pi * x) + rng.normal(scale=0.05, size=x.size)
    candidates = np.array([0.15, 0.25, 0.4])
    first = regression_bandwidth(
        x, y, candidates=candidates, n_folds=4, random_state=11
    )
    second = regression_bandwidth(
        x, y, candidates=candidates, n_folds=4, random_state=11
    )
    assert first in candidates
    assert first == second


def test_paired_bootstrap_band_is_reproducible_and_fixed_width():
    x = np.linspace(-1.0, 1.0, 35)
    y = np.sin(np.pi * x) + 0.1 * np.cos(7 * x)
    points = np.linspace(-0.8, 0.8, 11)
    first = regression_confidence_band(
        x, y, points, bandwidth=0.4, n_boot=20, random_state=19
    )
    second = regression_confidence_band(
        x, y, points, bandwidth=0.4, n_boot=20, random_state=19
    )
    np.testing.assert_array_equal(first.bootstrap_statistics, second.bootstrap_statistics)
    np.testing.assert_allclose(first.width, 2 * first.critical_value)


def test_mismatched_regression_inputs_are_rejected():
    with pytest.raises(ValueError, match="equal lengths"):
        debiased_local_linear([0, 1, 2, 3], [0, 1], bandwidth=1)
