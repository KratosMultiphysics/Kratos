---
type: C++ Method
title: RemoveSubModelPart
resource: sources/model_part.cpp#L2084-L2110
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
  - target: functions/sources/model_part/ModelPart/ErrorNonExistingSubModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPart::RemoveSubModelPart(std::string const& ThisSubModelPartName)`

# Calls

- [FullName](../../../../functions/includes/model_part/FullName.md)
- [GetSubModelPartNames](../../../../functions/sources/model_part/ModelPart/GetSubModelPartNames.md)
- [ErrorNonExistingSubModelPart](../../../../functions/sources/model_part/ModelPart/ErrorNonExistingSubModelPart.md)