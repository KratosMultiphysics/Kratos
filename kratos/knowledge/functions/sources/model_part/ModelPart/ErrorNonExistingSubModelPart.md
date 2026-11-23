---
type: C++ Method
title: ErrorNonExistingSubModelPart
resource: sources/model_part.cpp#L2397-L2407
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/GetSubModelPartNames
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part/ModelPart/pGetSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveSubModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPart::ErrorNonExistingSubModelPart(const std::string& rSubModelPartName) const`

# Calls

- [FullName](../../../../functions/includes/model_part/FullName.md)
- [GetSubModelPartNames](../../../../functions/sources/model_part/ModelPart/GetSubModelPartNames.md)

# Called by

- [pGetSubModelPart](../../../../functions/sources/model_part/ModelPart/pGetSubModelPart.md)
- [RemoveSubModelPart](../../../../functions/sources/model_part/ModelPart/RemoveSubModelPart.md)