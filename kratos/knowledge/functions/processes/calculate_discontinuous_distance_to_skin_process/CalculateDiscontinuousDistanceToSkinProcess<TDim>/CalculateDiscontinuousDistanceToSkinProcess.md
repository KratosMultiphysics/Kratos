---
type: C++ Method
title: CalculateDiscontinuousDistanceToSkinProcess
resource: processes/calculate_discontinuous_distance_to_skin_process.cpp#L38-L45
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`CalculateDiscontinuousDistanceToSkinProcess<TDim>::CalculateDiscontinuousDistanceToSkinProcess( ModelPart& rVolumePart, ModelPart& rSkinPart) : mFindIntersectedObjectsProcess(rVolumePart, rSkinPart) , mrSkinPart(rSkinPart) , mrVolumePart(rVolumePart)`

# Calls

- [RecursivelyValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)