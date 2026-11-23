---
type: C++ Method
title: AssignTablesToProperty
resource: utilities/read_materials_utility.cpp#L297-L333
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/CheckIfOverwritingTable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/properties/SetTable
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignMaterialToProperty
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ReadMaterialsUtility::AssignTablesToProperty( const Parameters MaterialData, Properties& rProperty )`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [CheckIfOverwritingTable](../../../../functions/utilities/read_materials_utility/CheckIfOverwritingTable.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [SetTable](../../../../functions/includes/properties/SetTable.md)

# Called by

- [AssignMaterialToProperty](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignMaterialToProperty.md)