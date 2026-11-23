---
type: C++ Method
title: ContraVariantPushForward
resource: sources/constitutive_law.cpp#L1356-L1365
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

`void ConstitutiveLaw::ContraVariantPushForward( Matrix& rMatrix, const Matrix& rF) //i.e. 2nd PK stress to Kirchhoff stress`

# Calls

- [prod](../../../../functions/includes/amatrix_interface/prod.md)