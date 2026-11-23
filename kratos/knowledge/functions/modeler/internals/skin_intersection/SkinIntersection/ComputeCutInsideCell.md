---
type: C++ Method
title: ComputeCutInsideCell
resource: modeler/internals/skin_intersection.h#L53-L71
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/coloring/color_and_intersect_with_cells_in_touch/ColorAndIntersectWithCellsInTouch/ApplyColorToCellsInTouchWithGeometryAndComputeIntersection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/entity_generation/generate_stl_intersection_with_cells/GenerateStlIntersectionWithCells/GenerateIntersection
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ComputeCutInsideCell( SplitGeometryPointerType pGeometry, const Point& rLowPoint, const Point& rHighPoint, std::size_t I, std::size_t J, std::size_t K)`

# Called by

- [ApplyColorToCellsInTouchWithGeometryAndComputeIntersection](../../../../../functions/modeler/coloring/color_and_intersect_with_cells_in_touch/ColorAndIntersectWithCellsInTouch/ApplyColorToCellsInTouchWithGeometryAndComputeIntersection.md)
- [GenerateIntersection](../../../../../functions/modeler/entity_generation/generate_stl_intersection_with_cells/GenerateStlIntersectionWithCells/GenerateIntersection.md)