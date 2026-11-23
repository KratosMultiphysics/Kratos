---
type: C++ Method
title: ConstitutiveMatrixTransformation
resource: sources/constitutive_law.cpp#L1426-L1474
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/sources/constitutive_law/ConstitutiveLaw/PullBackConstitutiveMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/constitutive_law/ConstitutiveLaw/PushForwardConstitutiveMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ConstitutiveLaw::ConstitutiveMatrixTransformation ( Matrix& rConstitutiveMatrix, const Matrix& rOriginalConstitutiveMatrix, const Matrix & rF )`

# Called by

- [PullBackConstitutiveMatrix](../../../../functions/sources/constitutive_law/ConstitutiveLaw/PullBackConstitutiveMatrix.md)
- [PushForwardConstitutiveMatrix](../../../../functions/sources/constitutive_law/ConstitutiveLaw/PushForwardConstitutiveMatrix.md)