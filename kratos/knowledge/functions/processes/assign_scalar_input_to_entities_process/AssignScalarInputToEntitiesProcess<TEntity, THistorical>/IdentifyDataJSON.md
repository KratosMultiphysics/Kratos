---
type: C++ Method
title: IdentifyDataJSON
resource: processes/assign_scalar_input_to_entities_process.cpp#L428-L476
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/IsVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/AssignScalarInputToEntitiesProcess
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AssignScalarInputToEntitiesProcess<TEntity, THistorical>::IdentifyDataJSON(const std::string& rFileName)`

# Calls

- [IsVector](../../../../functions/sources/kratos_parameters/Parameters/IsVector.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)

# Called by

- [AssignScalarInputToEntitiesProcess](../../../../functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/AssignScalarInputToEntitiesProcess.md)