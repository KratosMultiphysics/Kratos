---
type: C++ Method
title: ContraVariantPullBack
resource: sources/constitutive_law.cpp#L1372-L1384
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

`void ConstitutiveLaw::ContraVariantPullBack( Matrix& rMatrix, const Matrix& rF) //i.e. Kirchhoff stress to 2nd PK stress`

# Calls

- [prod](../../../../functions/includes/amatrix_interface/prod.md)