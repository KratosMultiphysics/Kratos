---
type: C++ Method
title: AddGeometry
resource: modeler/internals/cartesian_mesh_colors.h#L213-L245
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/internals/cartesian_mesh_colors/CalculateMinMaxNodePositions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CalculateMinMaxCenterOfElementPositions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_ray/CartesianRay/AddIntersection
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddGeometry(Element::GeometryType const& rGeometry, bool IsNodal)`

# Calls

- [CalculateMinMaxNodePositions](../../../../../functions/modeler/internals/cartesian_mesh_colors/CalculateMinMaxNodePositions.md)
- [CalculateMinMaxCenterOfElementPositions](../../../../../functions/modeler/internals/cartesian_mesh_colors/CalculateMinMaxCenterOfElementPositions.md)
- [AddIntersection](../../../../../functions/modeler/internals/cartesian_ray/CartesianRay/AddIntersection.md)