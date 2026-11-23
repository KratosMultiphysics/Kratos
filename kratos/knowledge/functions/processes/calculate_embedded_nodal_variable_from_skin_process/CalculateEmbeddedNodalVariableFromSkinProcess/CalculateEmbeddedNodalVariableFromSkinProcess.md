---
type: C++ Method
title: CalculateEmbeddedNodalVariableFromSkinProcess
resource: processes/calculate_embedded_nodal_variable_from_skin_process.h#L190-L198
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
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`CalculateEmbeddedNodalVariableFromSkinProcess( Model &rModel, Parameters rSettings) : CalculateEmbeddedNodalVariableFromSkinProcess( rModel.GetModelPart(rSettings["base_model_part_name"].GetString()), rModel.GetModelPart(rSettings["skin_model_part_name"].GetString()), [] (Parameters x) -> Parameters {x.ValidateAndAssignDefaults(StaticGetDefaultParameters()); return x;} (rSettings))`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)