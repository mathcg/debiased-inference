# Cross-language API contract

The R and Python packages use the same public function names and argument
defaults wherever the languages permit it:

| Capability | Function |
|---|---|
| Ordinary KDE bandwidth | `density_bandwidth` |
| Debiased KDE | `debiased_kde` |
| Density confidence band | `kde_confidence_band` |
| Density level set | `density_level_set` |
| Density level-set confidence set | `density_level_set_confidence` |
| Ordinary regression bandwidth | `regression_bandwidth` |
| Debiased local-linear regression | `debiased_local_linear` |
| Regression confidence band | `regression_confidence_band` |
| Inverse regression | `inverse_regression` |
| Inverse-regression confidence set | `inverse_regression_confidence` |
| Generic band inversion | `invert_confidence_band` |
| Finite-set Hausdorff distance | `hausdorff_distance` |

Both suites contain shared deterministic fixtures for the KDE and regression
estimators. These fixtures agree to at least 13 significant digits. Bootstrap
draws are reproducible within each language for a fixed `random_state`; exact
draw-by-draw equality across languages is not promised because R and NumPy use
different random-number generators.

`inverse_regression_confidence` accepts `method="normal"` in both languages
when the estimated inverse set has one crossing. It uses the bootstrap standard
deviation and the corresponding standard-normal critical value described in
Section 3.2.1 of the paper.
