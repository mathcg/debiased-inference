# debiased-inference (Python)

Python implementation of the procedures in Cheng and Chen,
“Nonparametric Inference via Bootstrapping the Debiased Estimator.”

```python
import numpy as np
from debiased_inference import kde_confidence_band

rng = np.random.default_rng(2026)
sample = rng.normal(size=300)
band = kde_confidence_band(sample, n_boot=499, random_state=2026)
print(band)
```

The package also provides debiased local-linear regression, density level-set
and inverse-regression confidence sets, normal-reference and cross-validated
bandwidth selectors, and a studentized density band. The project repository
contains the full statistical specification and cross-language API guide.
