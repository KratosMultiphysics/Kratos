---
type: C++ Method
title: ApplyConstraintsForPeriodicConditions
resource: processes/apply_periodic_boundary_condition_process.cpp#L157-L196
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/builtin_timer/BuiltinTimer/BuiltinTimer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/TransformNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/binbased_fast_point_locator/FindPointOnMeshSimplified
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyPeriodicConditionProcess::ApplyConstraintsForPeriodicConditions()`

# Calls

- [BuiltinTimer](../../../../functions/utilities/builtin_timer/BuiltinTimer/BuiltinTimer.md)
- [TransformNode](../../../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/TransformNode.md)
- [FindPointOnMeshSimplified](../../../../functions/utilities/binbased_fast_point_locator/FindPointOnMeshSimplified.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)