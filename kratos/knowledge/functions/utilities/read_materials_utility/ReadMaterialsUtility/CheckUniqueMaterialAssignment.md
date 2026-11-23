---
type: C++ Method
title: CheckUniqueMaterialAssignment
resource: utilities/read_materials_utility.cpp#L499-L543
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/CheckModelPartIsNotRepeated
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/GetPropertyBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ReadMaterialsUtility::CheckUniqueMaterialAssignment(Parameters Materials)`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [CheckModelPartIsNotRepeated](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/CheckModelPartIsNotRepeated.md)

# Called by

- [GetPropertyBlock](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/GetPropertyBlock.md)