---
type: C++ Function
title: HasInitialState
resource: includes/constitutive_law.h#L572-L575
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/constitutive_law/AddInitialStressVectorContribution
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/constitutive_law/AddInitialStrainVectorContribution
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/constitutive_law/AddInitialDeformationGradientMatrixContribution
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool HasInitialState() const`

# Called by

- [AddInitialStressVectorContribution](../../../functions/includes/constitutive_law/AddInitialStressVectorContribution.md)
- [AddInitialStrainVectorContribution](../../../functions/includes/constitutive_law/AddInitialStrainVectorContribution.md)
- [AddInitialDeformationGradientMatrixContribution](../../../functions/includes/constitutive_law/AddInitialDeformationGradientMatrixContribution.md)