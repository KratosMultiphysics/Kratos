---
type: C++ Method
title: ResetModelPart
resource: modeler/combine_model_part_modeler.cpp#L68-L77
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
  - target: functions/sources/model_part/ModelPart/RemoveElementsFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveConditionsFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CombineModelPartModeler::ResetModelPart(ModelPart& rCombinedModelPart) const`

# Calls

- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)
- [RemoveNodesFromAllLevels](../../../../functions/sources/model_part/ModelPart/RemoveNodesFromAllLevels.md)
- [RemoveElementsFromAllLevels](../../../../functions/sources/model_part/ModelPart/RemoveElementsFromAllLevels.md)
- [RemoveConditionsFromAllLevels](../../../../functions/sources/model_part/ModelPart/RemoveConditionsFromAllLevels.md)