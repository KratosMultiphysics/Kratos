---
type: C++ Method
title: AvoidZeroDistances
resource: processes/calculate_signed_distance_to_3d_skin_process.h#L1491-L1528
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcDistanceTo3DSkin
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AvoidZeroDistances( ModelPart::ElementsContainerType::iterator& Element, array_1d<double,4>& ElementalDistances)`

# Called by

- [CalcDistanceTo3DSkin](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcDistanceTo3DSkin.md)