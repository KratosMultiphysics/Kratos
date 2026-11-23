---
type: C++ Method
title: CoVariantPullBack
resource: sources/constitutive_law.cpp#L1410-L1420
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ConstitutiveLaw::CoVariantPullBack( Matrix& rMatrix, const Matrix& rF) //i.e. Almansi strain to Green-Lagrange strain`

# Calls

- [prod](../../../../functions/includes/amatrix_interface/prod.md)