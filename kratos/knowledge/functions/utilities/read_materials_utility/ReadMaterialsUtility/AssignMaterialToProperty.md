---
type: C++ Method
title: AssignMaterialToProperty
resource: utilities/read_materials_utility.cpp#L170-L190
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignVariablesToProperty
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignTablesToProperty
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignConstitutiveLawToProperty
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignAccessorsToProperty
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/CreateSubProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignPropertyBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ReadMaterialsUtility::AssignMaterialToProperty( const Parameters MaterialData, Properties& rProperty )`

# Calls

- [AssignVariablesToProperty](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignVariablesToProperty.md)
- [AssignTablesToProperty](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignTablesToProperty.md)
- [AssignConstitutiveLawToProperty](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignConstitutiveLawToProperty.md)
- [AssignAccessorsToProperty](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignAccessorsToProperty.md)

# Called by

- [CreateSubProperties](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/CreateSubProperties.md)
- [AssignPropertyBlock](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignPropertyBlock.md)