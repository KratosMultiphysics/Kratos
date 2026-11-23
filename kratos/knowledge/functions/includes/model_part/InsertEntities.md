---
type: C++ Function
title: InsertEntities
resource: includes/model_part.h#L2174-L2202
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/InsertEntityRange
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/includes/model_part/AddNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/AddMasterSlaveConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/AddElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/AddConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/AddGeometries
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void InsertEntities( const TContainerGetter& rContainerGetter, TInputContainer&& rInputContainer)`

# Calls

- [IsSubModelPart](../../../functions/includes/model_part/IsSubModelPart.md)
- [InsertEntityRange](../../../functions/includes/model_part/InsertEntityRange.md)

# Called by

- [AddNodes](../../../functions/includes/model_part/AddNodes.md)
- [AddMasterSlaveConstraints](../../../functions/includes/model_part/AddMasterSlaveConstraints.md)
- [AddElements](../../../functions/includes/model_part/AddElements.md)
- [AddConditions](../../../functions/includes/model_part/AddConditions.md)
- [AddGeometries](../../../functions/includes/model_part/AddGeometries.md)