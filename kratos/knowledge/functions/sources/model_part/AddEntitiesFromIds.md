---
type: C++ Function
title: AddEntitiesFromIds
resource: sources/model_part.cpp#L36-L84
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part/ModelPart/AddNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/AddElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/AddMasterSlaveConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/AddConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/AddGeometries
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddEntitiesFromIds( TContainerGetter&& rContainerGetter, ModelPart* pModelPart, const std::vector<IndexType>& rEntityIds)`

# Calls

- [IsSubModelPart](../../../functions/includes/model_part/IsSubModelPart.md)
- [FullName](../../../functions/includes/model_part/FullName.md)

# Called by

- [AddNodes](../../../functions/sources/model_part/ModelPart/AddNodes.md)
- [AddElements](../../../functions/sources/model_part/ModelPart/AddElements.md)
- [AddMasterSlaveConstraints](../../../functions/sources/model_part/ModelPart/AddMasterSlaveConstraints.md)
- [AddConditions](../../../functions/sources/model_part/ModelPart/AddConditions.md)
- [AddGeometries](../../../functions/sources/model_part/ModelPart/AddGeometries.md)