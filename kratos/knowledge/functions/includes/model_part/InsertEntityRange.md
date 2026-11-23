---
type: C++ Function
title: InsertEntityRange
resource: includes/model_part.h#L2129-L2171
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/IsSubSet
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
  - target: functions/includes/model_part/InsertEntities
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void InsertEntityRange( const TContainerGetter& rContainerGetter, TIterator begin, TIterator end)`

# Calls

- [IsSubModelPart](../../../functions/includes/model_part/IsSubModelPart.md)
- [IsSubSet](../../../functions/includes/model_part/IsSubSet.md)

# Called by

- [AddNodes](../../../functions/includes/model_part/AddNodes.md)
- [AddMasterSlaveConstraints](../../../functions/includes/model_part/AddMasterSlaveConstraints.md)
- [AddElements](../../../functions/includes/model_part/AddElements.md)
- [AddConditions](../../../functions/includes/model_part/AddConditions.md)
- [AddGeometries](../../../functions/includes/model_part/AddGeometries.md)
- [InsertEntities](../../../functions/includes/model_part/InsertEntities.md)