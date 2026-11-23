---
type: C++ Method
title: ApplyConstantVectorValueProcess
resource: processes/apply_constant_vectorvalue_process.cpp#L33-L39
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsDefined
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ApplyConstantVectorValueProcess::ApplyConstantVectorValueProcess( Model& rModel, Parameters ThisParameters ) : ApplyConstantVectorValueProcess(rModel.GetModelPart(ThisParameters["model_part_name"].GetString()), ThisParameters)`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [IsArray](../../../../functions/sources/kratos_parameters/Parameters/IsArray.md)
- [PrettyPrintJsonString](../../../../functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString.md)
- [IsNumber](../../../../functions/sources/kratos_parameters/Parameters/IsNumber.md)
- [IsString](../../../../functions/sources/kratos_parameters/Parameters/IsString.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)
- [IsDefined](../../../../functions/containers/flags/IsDefined.md)