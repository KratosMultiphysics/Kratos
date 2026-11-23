---
type: C++ Method
title: ApplyColorToCellsInTouchWithGeometryAndComputeIntersection
resource: modeler/coloring/color_and_intersect_with_cells_in_touch.cpp#L73-L97
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/internals/skin_intersection/ConvertToWorkGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CalculateOuterMinMaxNodePositions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/skin_intersection/SkinIntersection/ComputeCutInsideCell
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/coloring/color_and_intersect_with_cells_in_touch/ColorAndIntersectWithCellsInTouch/Apply
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ColorAndIntersectWithCellsInTouch::ApplyColorToCellsInTouchWithGeometryAndComputeIntersection( const Element::GeometryType& rGeometry, int InsideColor, Internals::SkinIntersection& rSkinIntersection) const`

# Calls

- [ConvertToWorkGeometry](../../../../../functions/modeler/internals/skin_intersection/ConvertToWorkGeometry.md)
- [CalculateOuterMinMaxNodePositions](../../../../../functions/modeler/internals/cartesian_mesh_colors/CalculateOuterMinMaxNodePositions.md)
- [GetElementalColor](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor.md)
- [ComputeCutInsideCell](../../../../../functions/modeler/internals/skin_intersection/SkinIntersection/ComputeCutInsideCell.md)

# Called by

- [Apply](../../../../../functions/modeler/coloring/color_and_intersect_with_cells_in_touch/ColorAndIntersectWithCellsInTouch/Apply.md)