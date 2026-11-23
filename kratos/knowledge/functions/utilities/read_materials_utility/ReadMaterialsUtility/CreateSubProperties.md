---
type: C++ Method
title: CreateSubProperties
resource: utilities/read_materials_utility.cpp#L364-L430
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RecursivelyHasProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignMaterialToProperty
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignPropertyBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ReadMaterialsUtility::CreateSubProperties( ModelPart& rModelPart, const Parameters SubPropertiesData, Properties& rProperty )`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [RecursivelyHasProperties](../../../../functions/sources/model_part/ModelPart/RecursivelyHasProperties.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [CreateNewProperties](../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)
- [AssignMaterialToProperty](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignMaterialToProperty.md)

# Called by

- [AssignPropertyBlock](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignPropertyBlock.md)