---
type: C++ Function
title: ResizePositions
resource: containers/variables_list.h#L484-L514
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/variables_list/GetHashIndex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/variable_data/SourceKey
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/containers/variables_list/SetPosition
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ResizePositions()`

# Calls

- [GetHashIndex](../../../functions/containers/variables_list/GetHashIndex.md)
- [SourceKey](../../../functions/containers/variable_data/SourceKey.md)

# Called by

- [SetPosition](../../../functions/containers/variables_list/SetPosition.md)