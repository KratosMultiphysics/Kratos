---
type: C++ Method
title: CalculateDistanceToSkinProcess
resource: processes/calculate_distance_to_skin_process.cpp#L28-L34
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`CalculateDistanceToSkinProcess<TDim>::CalculateDistanceToSkinProcess( ModelPart& rVolumePart, ModelPart& rSkinPart) : CalculateDiscontinuousDistanceToSkinProcess<TDim>(rVolumePart, rSkinPart)`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)