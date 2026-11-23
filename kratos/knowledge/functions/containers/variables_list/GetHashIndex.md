---
type: C++ Function
title: GetHashIndex
resource: containers/variables_list.h#L475-L477
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/containers/variables_list/Has
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/variables_list/SetPosition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/variables_list/GetPosition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/variables_list/ResizePositions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`SizeType GetHashIndex(std::size_t Key, std::size_t TableSize, std::size_t HashFunctionIndex) const`

# Called by

- [Has](../../../functions/containers/variables_list/Has.md)
- [SetPosition](../../../functions/containers/variables_list/SetPosition.md)
- [GetPosition](../../../functions/containers/variables_list/GetPosition.md)
- [ResizePositions](../../../functions/containers/variables_list/ResizePositions.md)