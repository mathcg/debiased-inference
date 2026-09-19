# API guide

The snippets below use Python syntax. R uses the same function and parameter
names, with `TRUE`/`FALSE` in place of Python booleans.

## Density estimation

```python
h = density_bandwidth(x, method="normal_reference")
h_cv = density_bandwidth(x, method="cv")
fit = debiased_kde(x, points=grid, bandwidth=h, tau=1.0)
band = kde_confidence_band(
    x, points=grid, bandwidth=h, tau=1.0,
    confidence=0.95, n_boot=999, random_state=2026,
)
variable_band = kde_confidence_band(
    x, points=grid, bandwidth=h, studentized=True,
    n_boot=999, random_state=2026,
)
```

For two-dimensional samples, pass an `(n, 2)` sample matrix and an `(m, 2)`
evaluation matrix. A confidence band is simultaneous across the supplied
points.

## Density level sets

```python
estimated_set = density_level_set(x, level=0.1, points=grid, bandwidth=h)
hausdorff_set = density_level_set_confidence(
    x, level=0.1, points=grid, bandwidth=h,
    method="hausdorff", n_boot=999, random_state=2026,
)
inverted_set = density_level_set_confidence(
    x, level=0.1, points=grid, bandwidth=h,
    method="inversion", n_boot=999, random_state=2026,
)
```

For 2-D Hausdorff inference, `grid` must contain every point of a rectangular
Cartesian product. The returned `geometry` (Python) or `roots` (R) is a contour
point matrix. `mask` identifies grid points in the confidence region.

## Regression and inverse regression

```python
h = regression_bandwidth(x, y, n_folds=5, random_state=2026)
fit = debiased_local_linear(x, y, points=grid, bandwidth=h)
band = regression_confidence_band(
    x, y, points=grid, bandwidth=h,
    confidence=0.95, n_boot=999, random_state=2026,
)
inverse = inverse_regression(x, y, level=target, points=grid, bandwidth=h)
inverse_set = inverse_regression_confidence(
    x, y, level=target, points=grid, bandwidth=h,
    method="hausdorff", n_boot=999, random_state=2026,
)
```

When the estimated inverse-regression set has exactly one crossing,
`method="normal"` implements the paper's normal interval with a bootstrap
standard-error estimate.

Use an interior regression grid when local cubic fits become unstable near
sparsely sampled boundaries, matching the paper's simulation practice.

## Results

Estimator results expose `points`, `estimate`, `bandwidth`, `tau`, and `method`.
Band results additionally expose `lower`, `upper`, `critical_value`,
`confidence`, `n_boot`, and the bootstrap statistics. Set results expose
`roots`/`geometry`, a grid `mask`, and (for Hausdorff inference) `radius`.
