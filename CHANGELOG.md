# Changelog

## 0.1.1 — 2026-09-19

- Validate density and regression confidence bands against the published EJS
  simulations in both R and Python.
- Match density defaults to the paper's Gaussian normal-scale bandwidth rule.
- Expand the default regression cross-validation range so high-curvature
  designs do not truncate the optimum.
- Clarify that the paper's tabulated band width is the bootstrap band radius.

## 0.1.0 — 2026-09-18

- Add Gaussian debiased KDE for one- and multi-dimensional samples.
- Add normal-reference and least-squares CV bandwidth selection.
- Add fixed-width and studentized empirical-bootstrap simultaneous density bands.
- Add debiased one-dimensional local-linear regression with K-fold bandwidth CV.
- Add paired-bootstrap simultaneous regression bands.
- Add one-dimensional roots, two-dimensional density contours, Hausdorff-bootstrap
  confidence sets, confidence-band inversion, and the paper's normal interval
  for a unique inverse-regression crossing.
- Add matching R and Python APIs, reproducible examples, and cross-language
  numerical fixtures.
