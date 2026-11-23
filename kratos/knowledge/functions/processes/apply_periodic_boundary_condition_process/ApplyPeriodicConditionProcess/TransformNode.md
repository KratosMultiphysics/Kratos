---
type: C++ Method
title: TransformNode
resource: processes/apply_periodic_boundary_condition_process.cpp#L288-L302
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/ApplyConstraintsForPeriodicConditions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyPeriodicConditionProcess::TransformNode(const array_1d<double, 3 >& rCoordinates, array_1d<double, 3 >& rTransformedCoordinates) const`

# Called by

- [ApplyConstraintsForPeriodicConditions](../../../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/ApplyConstraintsForPeriodicConditions.md)