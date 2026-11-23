---
type: C++ Method
title: ReadDataJSON
resource: processes/assign_scalar_input_to_entities_process.cpp#L565-L597
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/result_dabatase/SetCommonColumn
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/serializer/Serializer/GetVariableData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/AssignScalarInputToEntitiesProcess
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AssignScalarInputToEntitiesProcess<TEntity, THistorical>::ReadDataJSON(const std::string& rFileName)`

# Calls

- [SetCommonColumn](../../../../functions/utilities/result_dabatase/SetCommonColumn.md)
- [GetVariableData](../../../../functions/sources/serializer/Serializer/GetVariableData.md)

# Called by

- [AssignScalarInputToEntitiesProcess](../../../../functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/AssignScalarInputToEntitiesProcess.md)