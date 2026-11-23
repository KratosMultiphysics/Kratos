---
type: C++ Method
title: AssignPropertyBlock
resource: utilities/read_materials_utility.cpp#L435-L494
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RecursivelyHasProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/properties/HasVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/properties/HasTables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignMaterialToProperty
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/CreateSubProperties
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/GetPropertyBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ReadMaterialsUtility::AssignPropertyBlock(Parameters Data)`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [RecursivelyHasProperties](../../../../functions/sources/model_part/ModelPart/RecursivelyHasProperties.md)
- [HasVariables](../../../../functions/includes/properties/HasVariables.md)
- [HasTables](../../../../functions/includes/properties/HasTables.md)
- [CreateNewProperties](../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)
- [AssignMaterialToProperty](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignMaterialToProperty.md)
- [CreateSubProperties](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/CreateSubProperties.md)

# Called by

- [GetPropertyBlock](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/GetPropertyBlock.md)