---
type: C++ Function
title: SetIdNotGeneratedFromString
resource: geometries/geometry.h#L4224-L4227
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/geometry/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/GenerateSelfAssignedId
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline void SetIdNotGeneratedFromString(IndexType& Id)`

# Called by

- [Create](../../../functions/geometries/geometry/Create.md)
- [GenerateSelfAssignedId](../../../functions/geometries/geometry/GenerateSelfAssignedId.md)