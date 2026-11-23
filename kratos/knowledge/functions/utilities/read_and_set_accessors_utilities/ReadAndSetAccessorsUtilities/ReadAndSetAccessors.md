---
type: C++ Method
title: ReadAndSetAccessors
resource: utilities/read_and_set_accessors_utilities.cpp#L26-L60
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/kratos_components/pGet
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/properties/HasAccessor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/properties/SetAccessor
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignAccessorsToProperty
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ReadAndSetAccessorsUtilities::ReadAndSetAccessors( const Parameters MaterialData, Properties& rProperty )`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [pGet](../../../../functions/includes/kratos_components/pGet.md)
- [HasAccessor](../../../../functions/includes/properties/HasAccessor.md)
- [SetAccessor](../../../../functions/includes/properties/SetAccessor.md)

# Called by

- [AssignAccessorsToProperty](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignAccessorsToProperty.md)