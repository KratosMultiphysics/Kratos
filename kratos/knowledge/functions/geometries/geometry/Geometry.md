---
type: C++ Function
title: Geometry
resource: geometries/geometry.h#L227-L231
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/GenerateSelfAssignedId
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/GenerateId
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Geometry() : mId(GenerateSelfAssignedId()) , mpGeometryData(&GeometryDataInstance())`

# Calls

- [GenerateSelfAssignedId](../../../functions/geometries/geometry/GenerateSelfAssignedId.md)
- [GenerateId](../../../functions/geometries/geometry/GenerateId.md)