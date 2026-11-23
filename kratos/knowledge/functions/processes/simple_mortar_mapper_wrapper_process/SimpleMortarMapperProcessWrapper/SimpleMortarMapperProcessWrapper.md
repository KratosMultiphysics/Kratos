---
type: C++ Method
title: SimpleMortarMapperProcessWrapper
resource: processes/simple_mortar_mapper_wrapper_process.h#L92-L187
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`SimpleMortarMapperProcessWrapper( ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, Parameters ThisParameters = Parameters(R"({})" ), LinearSolverType::Pointer pThisLinearSolver = nullptr )`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [SetBool](../../../../functions/sources/kratos_parameters/Parameters/SetBool.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)