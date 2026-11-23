---
type: C++ Function
title: LocalOffset
resource: containers/variables_list_data_value_container.h#L866-L870
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/variable_data/SourceKey
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/containers/variables_list_data_value_container/VariablesListDataValueContainer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/variables_list_data_value_container/Resize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/variables_list_data_value_container/AssignZero
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/variables_list_data_value_container/DestructElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/variables_list_data_value_container/DestructAllElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/variables_list_data_value_container/AssignData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/variables_list_data_value_container/save
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/variables_list_data_value_container/load
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline SizeType LocalOffset(VariableData const & rThisVariable) const`

# Calls

- [SourceKey](../../../functions/containers/variable_data/SourceKey.md)

# Called by

- [VariablesListDataValueContainer](../../../functions/containers/variables_list_data_value_container/VariablesListDataValueContainer.md)
- [Resize](../../../functions/containers/variables_list_data_value_container/Resize.md)
- [AssignZero](../../../functions/containers/variables_list_data_value_container/AssignZero.md)
- [DestructElements](../../../functions/containers/variables_list_data_value_container/DestructElements.md)
- [DestructAllElements](../../../functions/containers/variables_list_data_value_container/DestructAllElements.md)
- [AssignData](../../../functions/containers/variables_list_data_value_container/AssignData.md)
- [save](../../../functions/containers/variables_list_data_value_container/save.md)
- [load](../../../functions/containers/variables_list_data_value_container/load.md)