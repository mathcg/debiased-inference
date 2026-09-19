"""Bootstrap inference with debiased nonparametric estimators."""

from ._results import ConfidenceBandResult, EstimateResult, SetEstimateResult
from .bandwidth import density_bandwidth, regression_bandwidth
from .kde import debiased_kde, kde_confidence_band
from .regression import debiased_local_linear, regression_confidence_band
from .sets import (
    density_level_set,
    density_level_set_confidence,
    hausdorff_distance,
    inverse_regression,
    inverse_regression_confidence,
    invert_confidence_band,
    level_set,
)

__all__ = [
    "ConfidenceBandResult",
    "EstimateResult",
    "SetEstimateResult",
    "debiased_kde",
    "debiased_local_linear",
    "density_bandwidth",
    "density_level_set",
    "density_level_set_confidence",
    "hausdorff_distance",
    "inverse_regression",
    "inverse_regression_confidence",
    "invert_confidence_band",
    "kde_confidence_band",
    "level_set",
    "regression_bandwidth",
    "regression_confidence_band",
]

__version__ = "0.1.1"
