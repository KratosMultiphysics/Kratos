---
type: C++ Method
title: Execute
resource: processes/entity_erase_process.cpp#L72-L87
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/variable_utils/SetFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveNodesFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EntitiesEraseProcess<Node>::Execute()`

# Calls

- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)
- [RemoveNodesFromAllLevels](../../../../functions/sources/model_part/ModelPart/RemoveNodesFromAllLevels.md)
- [RemoveNodes](../../../../functions/sources/model_part/ModelPart/RemoveNodes.md)