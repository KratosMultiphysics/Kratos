---
type: C++ Method
title: ApplyRayCastingProcess
resource: processes/apply_ray_casting_process.cpp#L28-L40
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
  - target: functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/CreateDistanceGetterFunctor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ApplyRayCastingProcess<TDim>::ApplyRayCastingProcess( ModelPart& rVolumePart, ModelPart& rSkinPart, Parameters ThisParameters) : mSettings(ThisParameters), mpFindIntersectedObjectsProcess(new FindIntersectedGeometricalObjectsProcess(rVolumePart, rSkinPart)), mIsSearchStructureAllocated(true)`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [CreateDistanceGetterFunctor](../../../../functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/CreateDistanceGetterFunctor.md)
- [SetDouble](../../../../functions/sources/kratos_parameters/Parameters/SetDouble.md)
- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)