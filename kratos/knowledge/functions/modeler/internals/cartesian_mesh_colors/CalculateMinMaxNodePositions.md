---
type: C++ Function
title: CalculateMinMaxNodePositions
resource: modeler/internals/cartesian_mesh_colors.h#L631-L644
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/internals/cartesian_mesh_colors/CalculateMinMaxPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateNodePosition
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/AddGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateMinMaxNodePositions(TPointsContainerType const& Points, array_1d< std::size_t, 3 >& MinNodePosition, array_1d< std::size_t, 3 >& MaxNodePosition)`

# Calls

- [CalculateMinMaxPoints](../../../../functions/modeler/internals/cartesian_mesh_colors/CalculateMinMaxPoints.md)
- [CalculateNodePosition](../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateNodePosition.md)

# Called by

- [AddGeometry](../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/AddGeometry.md)