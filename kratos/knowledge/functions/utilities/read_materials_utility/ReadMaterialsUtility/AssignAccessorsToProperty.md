---
type: C++ Method
title: AssignAccessorsToProperty
resource: utilities/read_materials_utility.cpp#L338-L348
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/read_and_set_accessors_utilities/ReadAndSetAccessorsUtilities/ReadAndSetAccessors
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignMaterialToProperty
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ReadMaterialsUtility::AssignAccessorsToProperty( const Parameters MaterialData, Properties& rProperty )`

# Calls

- [ReadAndSetAccessors](../../../../functions/utilities/read_and_set_accessors_utilities/ReadAndSetAccessorsUtilities/ReadAndSetAccessors.md)

# Called by

- [AssignMaterialToProperty](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignMaterialToProperty.md)