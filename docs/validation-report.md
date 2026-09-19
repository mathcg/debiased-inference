# Validation against the published simulations

This report compares the package with the published version of Cheng and Chen
(2019), [*Nonparametric Inference via Bootstrapping the Debiased
Estimator*](https://projecteuclid.org/journalArticle/Download?urlId=10.1214%2F19-EJS1575),
Electronic Journal of Statistics 13(1), doi:10.1214/19-EJS1575.

## Summary

The corrected R and Python implementations reproduce the paper's principal
density and regression findings. The numerical band radii closely match the
published tables, empirical coverage is compatible with the reported coverage
under Monte Carlo uncertainty, and the main-text experiments reproduce the
claimed contrasts between ordinary and debiased bands.

The validation found and corrected two default-selection problems that were
not detected by unit tests alone:

1. Density bandwidth selection used a robust `0.9 * min(sd, IQR / 1.34)` rule
   instead of the `ks` Gaussian normal-scale rule used in the paper. The package
   now uses `(4 / (3 n))^(1/5) * sd` in one dimension.
2. The default regression CV grid ended at its lower boundary in the paper's
   high-curvature design. Its range is now wide enough to contain the empirical
   optimum.

The example output also now calls the bootstrap critical value a *band radius*;
the full distance from the lower to upper endpoint is twice that quantity.

## Density benchmark: Appendix D, Model 2

The design is an equal mixture of `N(-1, (2/3)^2)` and `N(1, (2/3)^2)`, with
simultaneous coverage evaluated over `[-2, 2]`. The Python replication uses the
paper's full 1,000 Monte Carlo repetitions and 1,000 bootstrap repetitions.

| n | EJS coverage | Python coverage (95% MC CI) | R coverage (95% MC CI) | EJS radius | Python radius | R radius |
|---:|---:|---:|---:|---:|---:|---:|
| 500 | 0.956 | 0.948 (0.932–0.960) | 0.955 (0.917–0.976) | 0.066 | 0.0662 | 0.0659 |
| 1000 | 0.949 | 0.945 (0.929–0.958) | 0.950 (0.910–0.973) | 0.052 | 0.0518 | 0.0516 |
| 2000 | 0.953 | 0.947 (0.931–0.959) | 0.965 (0.930–0.983) | 0.040 | 0.0403 | 0.0402 |

The R run uses 200 Monte Carlo repetitions and 499 bootstrap repetitions. All
published coverages lie within the corresponding replication confidence
intervals, and radii agree to the displayed precision.

## Main-text density benchmark: Figure 4

The design is `0.6 N(0, 1) + 0.4 N(4, 1)`. Results below use 400 Monte Carlo
repetitions and 499 bootstrap repetitions at 95% confidence.

| n | bandwidth | ordinary coverage | debiased coverage | ordinary radius | debiased radius |
|---:|:---:|---:|---:|---:|---:|
| 500 | `h_RT / 2` | 0.905 | 0.970 | 0.0440 | 0.0629 |
| 500 | `h_RT` | 0.008 | 0.893 | 0.0249 | 0.0369 |
| 500 | `2 h_RT` | 0.000 | 0.000 | 0.0124 | 0.0201 |
| 1000 | `h_RT / 2` | 0.925 | 0.960 | 0.0344 | 0.0488 |
| 1000 | `h_RT` | 0.003 | 0.915 | 0.0199 | 0.0291 |
| 1000 | `2 h_RT` | 0.000 | 0.000 | 0.0103 | 0.0161 |
| 2000 | `h_RT / 2` | 0.895 | 0.940 | 0.0269 | 0.0379 |
| 2000 | `h_RT` | 0.003 | 0.940 | 0.0158 | 0.0230 |
| 2000 | `2 h_RT` | 0.000 | 0.000 | 0.0084 | 0.0129 |

This reproduces Figure 4's substantive conclusions: the ordinary bootstrap
severely undercovers at conventional and doubled bandwidths, while debiasing
restores coverage at the conventional bandwidth as sample size increases.
Undersmoothing works better, and oversmoothing remains invalid. It also
reproduces Figure 5's width ordering: at each sample size the debiased band at
`h_RT` is narrower than the ordinary band at `h_RT / 2`.

## Regression benchmark: Tables 1 and 2

For the paper's asymmetric regression function on `[-1, 1]`, the Python run
uses 100 Monte Carlo repetitions and 199 bootstrap repetitions. The R run is an
independent replication of the `n = 500` configuration.

| n | EJS coverage | Python coverage (95% MC CI) | EJS radius | Python radius |
|---:|---:|---:|---:|---:|
| 500 | 0.976 | 0.960 (0.902–0.984) | 0.090 | 0.0853 |
| 1000 | 0.976 | 0.990 (0.946–0.998) | 0.066 | 0.0622 |
| 2000 | 0.963 | 0.970 (0.915–0.990) | 0.049 | 0.0462 |

At `n = 500`, R gives coverage 0.950 (95% Monte Carlo interval 0.888–0.978)
and radius 0.0869. The published coverage lies inside both languages' Monte
Carlo intervals; the simulated radii are within about 5% of the published
values. A separate Figure 7 sine-regression run gives debiased coverage
0.970/0.920/0.950 at `h_CV / 2`, `h_CV`, and `2 h_CV`; every 95% Monte Carlo
interval contains nominal 0.95 coverage, reproducing the paper's robustness
claim.

## Additional correctness evidence

- R and Python deterministic KDE fixtures agree to at least 14 significant
  digits, including a direct closed-form Gaussian-kernel calculation.
- Numerical integration verifies that the effective debiased kernel has unit
  mass and zero second moment.
- R and Python local-polynomial fixtures agree to at least 13 significant
  digits, and cubic-polynomial second derivatives are recovered exactly to
  numerical tolerance.
- Package test suites, lint, static typing, and CRAN checks run independently
  in CI.

## Reproducibility and interpretation

Raw summaries are stored in `validation/results/`. Scripts, seeds, grid sizes,
and run commands are documented in `validation/README.md`. R and NumPy use
different random-number generators, so draw-by-draw agreement is neither
expected nor required. Coverage comparisons should use the reported Monte
Carlo intervals rather than exact equality to one historical simulation seed.
