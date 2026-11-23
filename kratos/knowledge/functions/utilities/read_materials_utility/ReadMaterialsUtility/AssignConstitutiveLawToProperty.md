---
type: C++ Method
title: AssignConstitutiveLawToProperty
resource: utilities/read_materials_utility.cpp#L195-L217
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignMaterialToProperty
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ReadMaterialsUtility::AssignConstitutiveLawToProperty( const Parameters MaterialData, Properties& rProperty )`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)

# Called by

- [AssignMaterialToProperty](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignMaterialToProperty.md)