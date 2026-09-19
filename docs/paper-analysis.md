# Paper analysis and implementation map

Source reviewed: arXiv:1702.07027v3 (63 pages), including the estimator
definitions, algorithms, assumptions, theorems, simulation study, astronomy
application, proofs, and supplementary simulations.

## Central idea

An ordinary second-order kernel smoother has bias of order `h^2`, which is not
negligible at its usual MISE/MSE-optimal bandwidth. Naively bootstrapping that
estimator therefore reproduces stochastic variation around a biased target and
can undercover. The paper subtracts an estimate of the leading bias using a
second-derivative estimate at `b = h / tau`. With fixed positive `tau`—usually
one—the residual bias becomes higher order while variance remains on the
ordinary estimator's order. Consequently, a conventional bandwidth selected
for the ordinary estimator acts like an inferentially valid undersmoothing
bandwidth for the debiased estimator.

The crucial practical point is that the derivative estimate need not be
consistent by itself. Debiasing trades leading bias for additional stochastic
variation, and the bootstrap captures that added variation.

## Implemented correspondence

| Paper component | Software component |
|---|---|
| Equation (3), debiased KDE | `debiased_kde` |
| Figure 2, fixed-width density band | `kde_confidence_band` |
| Remark 1, studentized density band | `kde_confidence_band(studentized=True/TRUE)` |
| Section 3.1.1, density level set | `density_level_set` |
| Hausdorff level-set confidence set | `density_level_set_confidence(method="hausdorff")` |
| Remark 2, band inversion | `density_level_set_confidence(method="inversion")` |
| Equation (5), debiased local linear | `debiased_local_linear` |
| Figure 3, paired-bootstrap regression band | `regression_confidence_band` |
| Section 3.2.1, inverse regression | `inverse_regression` |
| Inverse-regression confidence set | `inverse_regression_confidence` |

For a unique inverse-regression crossing, the latter also exposes the paper's
normal approximation with a bootstrap standard error through `method="normal"`.

## Assumptions users need to understand

The software computes finite-sample procedures; the paper's coverage guarantee
is asymptotic and depends on its regularity conditions.

- The kernel is symmetric, second order, sufficiently smooth, integrable, and
  belongs (with its derivatives) to an appropriate VC-type function class. The
  Gaussian kernel supplied here satisfies these kernel conditions.
- The density or regression function has `2 + delta` Hölder smoothness with
  `2/3 < delta <= 2`. The density is bounded and has compact support with the
  paper's boundary conditions.
- For density and regression bands, `tau` is fixed and positive. The bandwidth
  tends to zero, `n h^d / log(n)` tends to infinity for density (respectively
  `n h / log(n)` for regression), and `n h^(d+4) / log(n)` (respectively
  `n h^5 / log(n)`) has a finite nonnegative limit.
- Density level sets require the density gradient to be bounded away from zero
  on the target contour. Inverse regression analogously requires a nonzero
  regression derivative at target crossings.
- Regression assumes a one-dimensional covariate with positive continuous
  design density on a compact domain and a uniformly bounded conditional
  fourth moment of the response.

These conditions are not fully testable from one dataset. Package results must
therefore be interpreted as asymptotically justified estimates, not as a
finite-sample coverage certificate.

## Rates and interpretation

- Debiased KDE pointwise bias: `O(h^(2 + delta))`; variance: `O(1/(n h^d))`.
- Debiased local-linear pointwise bias at `h = O(n^(-1/5))`:
  `O(h^(2 + delta))`; variance: `O(1/(n h))`.
- The simultaneous bands use the bootstrap quantile of the maximum absolute
  error, not independent pointwise quantiles. This distinction is preserved in
  both packages.
- A fourth-order effective kernel may yield negative density estimates and
  negative lower confidence limits. The implementation deliberately does not
  clip them, because clipping changes the estimator and bootstrap statistic.

## Numerical interpretation

The paper defines suprema on continuous domains and exact equality level sets.
Software must approximate both on a finite evaluation grid. Finer grids give a
better approximation but increase computation. One-dimensional roots are
linearly interpolated between sign-changing adjacent grid values. Two-
dimensional density contours are interpolated on triangles of a complete
rectangular grid and represented as a finite point cloud for Hausdorff distance;
band inversion remains available as a stable alternative.

For regression, boundary regions can have ill-conditioned local cubic fits.
The examples follow the paper by restricting inference to an interior domain.
The implementations report a singular-fit error rather than silently returning
an invalid band. Bootstrap resamples with singular or empty-set estimates are
redrawn up to `max_attempts`.

## Source issues resolved explicitly

Section 3.2.1 displays the Hausdorff distance between the same inverse set on
both sides, which would be identically zero. The intended expression, supported
by the algorithm and theorem discussion, is the distance between the bootstrap
set and the original estimated set. The source also alternates between the
original and starred estimated set when writing the final Minkowski expansion.
The software consistently centers bootstrap distances on the original estimate
and expands the original estimate by the bootstrap quantile.

## Simulation findings relevant to defaults

The paper generally finds nominal simultaneous coverage for the debiased
procedure at conventional bandwidths, while naive bootstrap bands under-cover
unless the ordinary estimator is undersmoothed. It recommends cross-validation
over a simple rule of thumb for difficult density shapes. The packages
therefore provide least-squares CV for KDE and K-fold CV for regression while
retaining a fast density normal-reference default. Published reproduction
should pass the exact reported bandwidth explicitly.
