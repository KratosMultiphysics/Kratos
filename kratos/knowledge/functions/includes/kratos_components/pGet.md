---
type: C++ Function
title: pGet
resource: includes/kratos_components.h#L343-L348
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/kratos_components/GetMessageUnregisteredVariable
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/containers/data_value_container/DataValueContainer/load
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/periodic_variables_container/PeriodicVariablesContainer/load
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/variables_list/VariablesList/load
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/serializer/Serializer/GetVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/table_accessor/TableAccessor/load
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_and_set_accessors_utilities/ReadAndSetAccessorsUtilities/ReadAndSetAccessors
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static const VariableData* pGet(const std::string& rName)`

# Calls

- [GetMessageUnregisteredVariable](../../../functions/includes/kratos_components/GetMessageUnregisteredVariable.md)

# Called by

- [load](../../../functions/containers/data_value_container/DataValueContainer/load.md)
- [load](../../../functions/containers/periodic_variables_container/PeriodicVariablesContainer/load.md)
- [load](../../../functions/containers/variables_list/VariablesList/load.md)
- [GetVariableData](../../../functions/sources/serializer/Serializer/GetVariableData.md)
- [load](../../../functions/sources/table_accessor/TableAccessor/load.md)
- [ReadAndSetAccessors](../../../functions/utilities/read_and_set_accessors_utilities/ReadAndSetAccessorsUtilities/ReadAndSetAccessors.md)