# debiasedInference

R implementation of the procedures in Cheng and Chen,
[“Nonparametric Inference via Bootstrapping the Debiased
Estimator”](https://projecteuclid.org/journalArticle/Download?urlId=10.1214%2F19-EJS1575).

If you use these methods or this software, please cite the paper (Electronic
Journal of Statistics 13(1), 2019; doi:10.1214/19-EJS1575).

```r
library(debiasedInference)

set.seed(2026)
x <- rnorm(300)
fit <- debiased_kde(x)
band <- kde_confidence_band(x, n_boot = 499, random_state = 2026)
```

The result objects are ordinary lists with classes `di_estimate`, `di_band`,
and `di_set`, so their numerical components are directly accessible.
