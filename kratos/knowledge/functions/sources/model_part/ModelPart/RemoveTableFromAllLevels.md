---
type: C++ Method
title: RemoveTableFromAllLevels
resource: sources/model_part.cpp#L591-L600
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveTable
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPart::RemoveTableFromAllLevels(ModelPart::IndexType TableId)`

# Calls

- [IsSubModelPart](../../../../functions/includes/model_part/IsSubModelPart.md)
- [RemoveTable](../../../../functions/sources/model_part/ModelPart/RemoveTable.md)