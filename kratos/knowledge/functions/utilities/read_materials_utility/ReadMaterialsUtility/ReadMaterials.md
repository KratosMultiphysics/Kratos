---
type: C++ Method
title: ReadMaterials
resource: utilities/read_materials_utility.cpp#L107-L114
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/GetPropertyBlock
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_serializer/create_model_part
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ReadMaterialsUtility::ReadMaterials(Parameters MaterialData)`

# Calls

- [GetPropertyBlock](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/GetPropertyBlock.md)

# Called by

- [create_model_part](../../../../functions/tests/test_serializer/create_model_part.md)