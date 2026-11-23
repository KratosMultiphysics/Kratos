---
type: C++ Method
title: CoVariantPushForward
resource: sources/constitutive_law.cpp#L1391-L1403
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

`void ConstitutiveLaw::CoVariantPushForward( Matrix& rMatrix, const Matrix& rF) //i.e. Green-Lagrange strain to Almansi strain`

# Calls

- [prod](../../../../functions/includes/amatrix_interface/prod.md)