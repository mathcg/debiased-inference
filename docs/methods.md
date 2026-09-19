# Statistical specification

This document is the shared source of truth for the R and Python packages.
Notation follows Cheng and Chen, *Nonparametric Inference via Bootstrapping the
Debiased Estimator*, arXiv:1702.07027v3.

## Debiased kernel density estimation

For observations in `d` dimensions, scalar bandwidth `h`, and `tau = h / b`,
the estimator is

```text
p_hat_tau,h(z) = (n h^d)^(-1) sum_i M_tau((z - X_i) / h)
M_tau(u) = K(u) - 0.5 c_K tau^(d+2) Laplacian(K)(tau u).
```

For the standard multivariate Gaussian kernel, `c_K = 1` and
`Laplacian(K)(u) = (||u||^2 - d) K(u)`. The default `tau` is 1, as recommended
in the paper. The estimator may be negative because `M_tau` is a fourth-order
kernel; values are not clipped.

The fixed-width simultaneous band resamples observations with replacement,
recomputes the estimator using the original bandwidth and grid, and takes the
empirical `(1 - alpha)` quantile of
`max_grid |p_hat_star - p_hat|`. The studentized band uses the finite-sample
variance expression in Remark 1 of the paper. Quantiles use the standard
linearly interpolated empirical quantile in both languages.

## Debiased local-linear regression

At each evaluation point, the ordinary local-linear estimate is the intercept
of a weighted least-squares fit of degree one. The second derivative estimate
is `2!` times the quadratic coefficient of a weighted degree-three fit with
bandwidth `b = h / tau`. The debiased estimate is

```text
r_hat_tau,h(z) = r_hat_h(z) - 0.5 c_K h^2 r_hat_b^(2)(z).
```

Regression confidence bands use the paired empirical bootstrap and the
`L-infinity` distance on the evaluation grid, exactly as in Figure 3 of the
paper. Singular local fits produce missing estimates; bootstrap replicates with
any missing grid estimate are rejected and redrawn up to a documented limit.

## Bandwidths

Bandwidth selection is performed for the *ordinary* estimator, not the
debiased estimator. Density supports a normal-reference rule (the default) and
least-squares cross-validation. Regression defaults to deterministic K-fold
cross-validation of the ordinary local-linear smoother. Explicit positive
bandwidths are always accepted and are preferred when reproducing a published
analysis.

## Numerical contract

- Inputs must be finite and contain at least two distinct observations.
- A bandwidth is a finite positive scalar. The current release implements the
  isotropic bandwidth used in the paper; bandwidth matrices are out of scope.
- Evaluation grids are sorted only when generated automatically. User-supplied
  order is preserved.
- Randomness is local to the function call and controlled by `random_state`.
- Confidence levels are strictly between zero and one; `n_boot >= 1`.
- A returned band covers only the supplied finite grid. Continuous-domain
  coverage is approximated by making that grid sufficiently dense.
- Equality level sets support one-dimensional grids and complete rectangular
  two-dimensional grids. Two-dimensional contours are represented by a finite
  cloud of linearly interpolated contour points.

## Paper erratum applied in software

The inverse-regression bootstrap display in Section 3.2.1 writes the Hausdorff
distance between the same set twice. The operative algorithm and surrounding
text require the distance between the bootstrap set and the original estimate;
the packages implement that intended expression.
