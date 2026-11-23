---
type: C++ Method
title: AssignScalarFieldToEntitiesProcess
resource: processes/assign_scalar_field_to_entities_process.cpp#L26-L34
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
---

# Signature

`AssignScalarFieldToEntitiesProcess<TEntity, THistorical>::AssignScalarFieldToEntitiesProcess( Model& rModel, Parameters ThisParameters ) : AssignScalarFieldToEntitiesProcess(rModel.GetModelPart(ThisParameters["model_part_name"].GetString()), ThisParameters)`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)