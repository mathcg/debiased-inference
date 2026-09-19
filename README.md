# debiasedInference

[![CI](https://github.com/mathcg/debiased-inference/actions/workflows/ci.yml/badge.svg)](https://github.com/mathcg/debiased-inference/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/debiased-inference.svg)](https://pypi.org/project/debiased-inference/)

Reference R and Python implementations for Gang Cheng and Yen-Chi Chen,
[“Nonparametric Inference via Bootstrapping the Debiased
Estimator”](https://projecteuclid.org/journalArticle/Download?urlId=10.1214%2F19-EJS1575).

The project implements the paper's debiased kernel density estimator and
debiased local-linear regression estimator, together with empirical-bootstrap
simultaneous confidence bands. The R and Python packages intentionally expose
matching function names and defaults.

See the [paper analysis](docs/paper-analysis.md), [statistical
specification](docs/methods.md), and [API guide](docs/api.md) for the reasoning,
implementation conventions, and complete workflow.

## Layout

- `python/`: installable Python package
- `r/`: installable R package
- `docs/`: mathematical and implementation notes
- `examples/`: reproducible cross-language examples

The `paper_density_simulation` scripts reproduce the paper's two-component
Gaussian-mixture coverage design at configurable Monte Carlo and bootstrap
sizes; their headers show the paper-scale settings.

## Python quick start

```console
python -m pip install debiased-inference
```

```python
import numpy as np
from debiased_inference import kde_confidence_band

rng = np.random.default_rng(2026)
x = rng.normal(size=300)
band = kde_confidence_band(x, n_boot=499, random_state=2026)
```

## R quick start

```r
install.packages(
  "https://github.com/mathcg/debiased-inference/releases/download/v0.1.0/debiasedInference_0.1.0.tar.gz",
  repos = NULL,
  type = "source"
)
library(debiasedInference)

set.seed(2026)
x <- rnorm(300)
band <- kde_confidence_band(x, n_boot = 499, random_state = 2026)
```

## Scope

Version `0.1.0` implements the Gaussian-kernel procedures analyzed and used in the
paper. Confidence bands are simultaneous over the supplied evaluation grid;
the grid is therefore part of the numerical approximation and should cover the
scientific domain of interest densely.

## Citation

If you use these methods or this software in your work, please cite the paper:

> Cheng, G. and Chen, Y.-C. (2019). Nonparametric inference via bootstrapping
> the debiased estimator. *Electronic Journal of Statistics*, 13(1).
> https://doi.org/10.1214/19-EJS1575

The [published paper on Project
Euclid](https://projecteuclid.org/journalArticle/Download?urlId=10.1214%2F19-EJS1575)
is the canonical methodological reference for this repository. Citation
metadata is also available in [`CITATION.cff`](CITATION.cff).
