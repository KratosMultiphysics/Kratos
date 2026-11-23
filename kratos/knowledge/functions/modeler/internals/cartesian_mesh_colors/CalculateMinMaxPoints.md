---
type: C++ Function
title: CalculateMinMaxPoints
resource: modeler/internals/cartesian_mesh_colors.h#L647-L660
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/internals/cartesian_mesh_colors/CalculateOuterMinMaxNodePositions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CalculateMinMaxNodePositions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateMinMaxPoints(TPointsContainerType const& Points, Point& MinPoint, Point& MaxPoint) const`

# Called by

- [CalculateOuterMinMaxNodePositions](../../../../functions/modeler/internals/cartesian_mesh_colors/CalculateOuterMinMaxNodePositions.md)
- [CalculateMinMaxNodePositions](../../../../functions/modeler/internals/cartesian_mesh_colors/CalculateMinMaxNodePositions.md)