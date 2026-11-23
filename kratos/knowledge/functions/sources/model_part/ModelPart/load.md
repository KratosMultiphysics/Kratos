---
type: C++ Method
title: load
resource: sources/model_part.cpp#L2349-L2394
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/serializer/IsDataOnly
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/SetParentModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPart::load(Serializer& rSerializer)`

# Calls

- [IsSubModelPart](../../../../functions/includes/model_part/IsSubModelPart.md)
- [IsDataOnly](../../../../functions/includes/serializer/IsDataOnly.md)
- [SetParentModelPart](../../../../functions/includes/model_part/SetParentModelPart.md)