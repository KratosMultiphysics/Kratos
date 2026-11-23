---
type: C++ Method
title: AssignScalarVariableToEntitiesProcess
resource: processes/assign_scalar_variable_to_entities_process.cpp#L27-L35
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
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`AssignScalarVariableToEntitiesProcess<TEntity, THistorical>::AssignScalarVariableToEntitiesProcess( Model& rModel, Parameters ThisParameters ) : AssignScalarVariableToEntitiesProcess(rModel.GetModelPart(ThisParameters["model_part_name"].GetString()), ThisParameters)`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)