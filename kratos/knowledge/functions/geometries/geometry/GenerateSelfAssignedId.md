---
type: C++ Function
title: GenerateSelfAssignedId
resource: geometries/geometry.h#L4197-L4209
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/SetIdSelfAssigned
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/SetIdNotGeneratedFromString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/geometry/Geometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`IndexType GenerateSelfAssignedId() const`

# Calls

- [SetIdSelfAssigned](../../../functions/geometries/geometry/SetIdSelfAssigned.md)
- [SetIdNotGeneratedFromString](../../../functions/geometries/geometry/SetIdNotGeneratedFromString.md)

# Called by

- [Geometry](../../../functions/geometries/geometry/Geometry.md)