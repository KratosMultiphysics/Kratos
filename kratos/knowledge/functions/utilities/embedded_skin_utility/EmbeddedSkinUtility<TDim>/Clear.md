---
type: C++ Method
title: Clear
resource: utilities/embedded_skin_utility.cpp#L101-L115
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/RemoveNodesFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveElementsFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveConditionsFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EmbeddedSkinUtility<TDim>::Clear()`

# Calls

- [RemoveNodesFromAllLevels](../../../../functions/sources/model_part/ModelPart/RemoveNodesFromAllLevels.md)
- [RemoveElementsFromAllLevels](../../../../functions/sources/model_part/ModelPart/RemoveElementsFromAllLevels.md)
- [RemoveConditionsFromAllLevels](../../../../functions/sources/model_part/ModelPart/RemoveConditionsFromAllLevels.md)