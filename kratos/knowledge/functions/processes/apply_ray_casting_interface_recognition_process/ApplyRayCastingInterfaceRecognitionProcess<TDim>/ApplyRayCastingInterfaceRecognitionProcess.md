---
type: C++ Method
title: ApplyRayCastingInterfaceRecognitionProcess
resource: processes/apply_ray_casting_interface_recognition_process.cpp#L18-L31
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/CreateDistanceGetterFunctor
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ApplyRayCastingInterfaceRecognitionProcess<TDim>::ApplyRayCastingInterfaceRecognitionProcess( Model& rModel, Parameters ThisParameters) : ApplyRayCastingProcess<TDim>( rModel.GetModelPart(ThisParameters["volume_model_part"].GetString()), rModel.GetModelPart(ThisParameters["skin_model_part"].GetString()), BaseType::GetDefaultParameters())`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [CreateDistanceGetterFunctor](../../../../functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/CreateDistanceGetterFunctor.md)