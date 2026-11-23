---
type: C++ Function
title: IsSlip
resource: utilities/coordinate_transformation_utilities.h#L1008-L1011
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/coordinate_transformation_utilities/ApplySlipCondition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/coordinate_transformation_utilities/RotateVelocities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/coordinate_transformation_utilities/RecoverVelocities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/coordinate_transformation_utilities/RotateAuxPure
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/coordinate_transformation_utilities/RotateRHSAux
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool IsSlip(const Node& rNode) const`

# Called by

- [ApplySlipCondition](../../../functions/utilities/coordinate_transformation_utilities/ApplySlipCondition.md)
- [RotateVelocities](../../../functions/utilities/coordinate_transformation_utilities/RotateVelocities.md)
- [RecoverVelocities](../../../functions/utilities/coordinate_transformation_utilities/RecoverVelocities.md)
- [RotateAuxPure](../../../functions/utilities/coordinate_transformation_utilities/RotateAuxPure.md)
- [RotateRHSAux](../../../functions/utilities/coordinate_transformation_utilities/RotateRHSAux.md)