---
type: C++ Method
title: Execute
resource: processes/entity_erase_process.cpp#L93-L108
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/variable_utils/SetFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveElementsFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveElements
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EntitiesEraseProcess<Element>::Execute()`

# Calls

- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)
- [RemoveElementsFromAllLevels](../../../../functions/sources/model_part/ModelPart/RemoveElementsFromAllLevels.md)
- [RemoveElements](../../../../functions/sources/model_part/ModelPart/RemoveElements.md)