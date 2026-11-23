---
type: C++ Method
title: Update
resource: solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip.h#L170-L185
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/coordinate_transformation_utilities/RotateVelocities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/coordinate_transformation_utilities/RecoverVelocities
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Update(ModelPart& r_model_part, DofsArrayType& rDofSet, TSystemMatrixType& A, TSystemVectorType& Dx, TSystemVectorType& b) override`

# Calls

- [RotateVelocities](../../../../../functions/utilities/coordinate_transformation_utilities/RotateVelocities.md)
- [RecoverVelocities](../../../../../functions/utilities/coordinate_transformation_utilities/RecoverVelocities.md)