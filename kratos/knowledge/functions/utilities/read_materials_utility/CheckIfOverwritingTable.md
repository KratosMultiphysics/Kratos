---
type: C++ Function
title: CheckIfOverwritingTable
resource: utilities/read_materials_utility.cpp#L41-L49
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/properties/HasTable
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignTablesToProperty
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CheckIfOverwritingTable(const Properties& rProps, const Variable<TValueType>& rInputVariable, const Variable<TValueType>& rOutputVariable)`

# Calls

- [HasTable](../../../functions/includes/properties/HasTable.md)

# Called by

- [AssignTablesToProperty](../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignTablesToProperty.md)