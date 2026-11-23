---
type: C++ Method
title: FindGlobalNodalEntityNeighboursProcess
resource: processes/find_global_nodal_entity_neighbours_process.cpp#L72-L80
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/project/Project/GetModel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`FindGlobalNodalEntityNeighboursProcess<TContainerType>::FindGlobalNodalEntityNeighboursProcess( Model& rModel, Parameters Params) : mrModel(rModel), mrOutputVariable(GetDefaultOutputVariable())`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetModel](../../../../functions/python_scripts/project/Project/GetModel.md)
- [FullName](../../../../functions/includes/model_part/FullName.md)