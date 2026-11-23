---
type: C++ Method
title: ResetModelPart
resource: modeler/connectivity_preserve_modeler.cpp#L207-L222
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
  - target: functions/sources/model_part/ModelPart/RemoveMasterSlaveConstraintsFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ConnectivityPreserveModeler::ResetModelPart(ModelPart& rDestinationModelPart) const`

# Calls

- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)
- [RemoveNodesFromAllLevels](../../../../functions/sources/model_part/ModelPart/RemoveNodesFromAllLevels.md)
- [RemoveElementsFromAllLevels](../../../../functions/sources/model_part/ModelPart/RemoveElementsFromAllLevels.md)
- [RemoveConditionsFromAllLevels](../../../../functions/sources/model_part/ModelPart/RemoveConditionsFromAllLevels.md)
- [RemoveMasterSlaveConstraintsFromAllLevels](../../../../functions/sources/model_part/ModelPart/RemoveMasterSlaveConstraintsFromAllLevels.md)