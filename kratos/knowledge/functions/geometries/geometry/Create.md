---
type: C++ Function
title: Create
resource: geometries/geometry.h#L875-L895
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
  - target: functions/geometries/geometry/SetIdWithoutCheck
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/gid_io/GidIOBase/GetData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual Pointer Create( PointsArrayType const& rThisPoints ) const`

# Calls

- [SetIdSelfAssigned](../../../functions/geometries/geometry/SetIdSelfAssigned.md)
- [SetIdNotGeneratedFromString](../../../functions/geometries/geometry/SetIdNotGeneratedFromString.md)
- [SetIdWithoutCheck](../../../functions/geometries/geometry/SetIdWithoutCheck.md)
- [GetData](../../../functions/sources/gid_io/GidIOBase/GetData.md)