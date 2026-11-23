---
type: C++ Method
title: GetPropertyBlock
resource: utilities/read_materials_utility.cpp#L119-L150
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/CheckUniqueMaterialAssignment
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RecursivelyHasProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignPropertyBlock
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/ReadMaterialsUtility
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/ReadMaterials
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ReadMaterialsUtility::GetPropertyBlock(Parameters Materials)`

# Calls

- [CheckUniqueMaterialAssignment](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/CheckUniqueMaterialAssignment.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [RecursivelyHasProperties](../../../../functions/sources/model_part/ModelPart/RecursivelyHasProperties.md)
- [CreateNewProperties](../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)
- [AssignPropertyBlock](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignPropertyBlock.md)

# Called by

- [ReadMaterialsUtility](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/ReadMaterialsUtility.md)
- [ReadMaterials](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/ReadMaterials.md)