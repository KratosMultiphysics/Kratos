---
type: C++ Method
title: CalculateElementalRayColors
resource: modeler/internals/cartesian_mesh_colors.h#L309-L386
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/internals/cartesian_ray/CartesianRay/CollapseIntersectionPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_ray/CartesianRay/CalculateColor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/coloring/color_cells_with_inside_center/ColorCellsWithInsideCenter/Apply
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateElementalRayColors(array_1d< std::size_t, 3 > const& MinRayPosition, array_1d< std::size_t, 3 > const& MaxRayPosition, int InsideColor, int OutsideColor)`

# Calls

- [CollapseIntersectionPoints](../../../../../functions/modeler/internals/cartesian_ray/CartesianRay/CollapseIntersectionPoints.md)
- [CalculateColor](../../../../../functions/modeler/internals/cartesian_ray/CartesianRay/CalculateColor.md)
- [GetElementalColor](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor.md)

# Called by

- [Apply](../../../../../functions/modeler/coloring/color_cells_with_inside_center/ColorCellsWithInsideCenter/Apply.md)