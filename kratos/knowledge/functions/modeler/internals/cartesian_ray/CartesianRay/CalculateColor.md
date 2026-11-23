---
type: C++ Method
title: CalculateColor
resource: modeler/internals/cartesian_ray.h#L138-L169
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateNodalRayColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateElementalRayColors
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateColor(std::vector<double> const& Coordinates, int InsideColor, int OutsideColor, std::vector<double>& ResultingColors, double NearEnough)`

# Called by

- [CalculateNodalRayColors](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateNodalRayColors.md)
- [CalculateElementalRayColors](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateElementalRayColors.md)