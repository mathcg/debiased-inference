# Contributing

Changes to statistical formulas should update `docs/methods.md`, both language
implementations, and their shared numerical fixtures.

## Python checks

```bash
cd python
uv run --extra test pytest --cov=debiased_inference
uv run --extra dev ruff check src tests
uv build
```

## R checks

From the repository root with R and `testthat` installed:

```bash
R CMD build r
R CMD check --no-manual debiasedInference_*.tar.gz
```

Do not change an estimator to clip negative density values: negative values are
an expected consequence of the fourth-order effective kernel. New stochastic
tests must set `random_state`; deterministic estimator tests should use explicit
bandwidths.
