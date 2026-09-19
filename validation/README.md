# Published-simulation validation

This directory validates the package against the simulations in Cheng and
Chen (2019), *Nonparametric Inference via Bootstrapping the Debiased
Estimator*, Electronic Journal of Statistics 13(1),
doi:10.1214/19-EJS1575.

`paper_density_validation.py` reproduces two published designs using the public
Python package API:

- `--study model2` corresponds to Appendix D, Tables 3 and 4: the equal mixture
  of `N(-1, (2/3)^2)` and `N(1, (2/3)^2)`, evaluated on `[-2, 2]`.
- `--study main` corresponds to Figure 4: the mixture `0.6 N(0, 1) +
  0.4 N(4, 1)`. It can compare the package's debiased band with an ordinary-KDE
  bootstrap reference at `h_RT / 2`, `h_RT`, and `2 h_RT`.

The paper used 1,000 Monte Carlo repetitions and 1,000 bootstrap repetitions.
Run its scale with:

```console
uv run --project python python validation/paper_density_validation.py \
  --study model2 --replications 1000 --n-boot 1000 --workers 4 \
  --output validation/results/model2-python.csv
```

The reported `mean_band_radius` is the bootstrap critical value, matching the
quantity called average band width in Tables 4 and 6. The geometric full width
between the lower and upper endpoints is twice this value.

`paper_density_validation.R` independently runs the Appendix D Model 2 design
through the public R API. The manually dispatched `Validate published
simulations` GitHub workflow runs it on a clean Linux R installation and saves
the summary as a workflow artifact.

`paper_regression_validation.py` reproduces the sine-regression design from
Figure 7 or the published debiased-CV column from Tables 1 and 2. It uses the
public `regression_bandwidth` and `regression_confidence_band` APIs throughout.
`paper_regression_validation.R` provides an independent R replication of the
`n = 500` Table 1 configuration.
