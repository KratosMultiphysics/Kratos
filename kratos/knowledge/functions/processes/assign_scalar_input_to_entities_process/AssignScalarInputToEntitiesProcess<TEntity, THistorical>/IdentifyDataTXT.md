---
type: C++ Method
title: IdentifyDataTXT
resource: processes/assign_scalar_input_to_entities_process.cpp#L368-L422
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/string_utilities/ContainsPartialString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/string_utilities/ErasePartialString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/AssignScalarInputToEntitiesProcess
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AssignScalarInputToEntitiesProcess<TEntity, THistorical>::IdentifyDataTXT(const std::string& rFileName)`

# Calls

- [ContainsPartialString](../../../../functions/utilities/string_utilities/ContainsPartialString.md)
- [ErasePartialString](../../../../functions/utilities/string_utilities/ErasePartialString.md)

# Called by

- [AssignScalarInputToEntitiesProcess](../../../../functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/AssignScalarInputToEntitiesProcess.md)