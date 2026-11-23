---
type: C++ Method
title: AssignVariablesToProperty
resource: utilities/read_materials_utility.cpp#L222-L292
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/FilterVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/CheckIfOverwritingValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignMaterialToProperty
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ReadMaterialsUtility::AssignVariablesToProperty( const Parameters MaterialData, Properties& rProperty )`

# Calls

- [FilterVariables](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/FilterVariables.md)
- [CheckIfOverwritingValue](../../../../functions/utilities/read_materials_utility/CheckIfOverwritingValue.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)

# Called by

- [AssignMaterialToProperty](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignMaterialToProperty.md)