---
type: Python Function
title: create_model_part
resource: tests/test_serializer.py#L108-L170
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/ReadMaterials
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_serializer/TestSerializer/test_DataOnly
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def create_model_part(factor = 1) -> KratosMultiphysics.Model:`

# Calls

- [ReadMaterials](../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/ReadMaterials.md)

# Called by

- [test_DataOnly](../../../functions/tests/test_serializer/TestSerializer/test_DataOnly.md)