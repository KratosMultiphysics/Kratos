---
type: C++ Function
title: AddInitialDeformationGradientMatrixContribution
resource: includes/constitutive_law.h#L605-L611
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/constitutive_law/HasInitialState
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddInitialDeformationGradientMatrixContribution(TMatrixType& rF)`

# Calls

- [HasInitialState](../../../functions/includes/constitutive_law/HasInitialState.md)
- [prod](../../../functions/includes/amatrix_interface/prod.md)