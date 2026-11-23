---
type: C++ Method
title: ApplyColorToCellsInTouchWithGeometry
resource: modeler/coloring/color_cells_in_touch.cpp#L60-L81
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/internals/cartesian_mesh_colors/CalculateOuterMinMaxNodePositions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/coloring/color_cells_in_touch/ColorCellsInTouch/Apply
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ColorCellsInTouch::ApplyColorToCellsInTouchWithGeometry( const Geometry<Node>& rGeometry, int InsideColor, int OutsideColor) const`

# Calls

- [CalculateOuterMinMaxNodePositions](../../../../../functions/modeler/internals/cartesian_mesh_colors/CalculateOuterMinMaxNodePositions.md)
- [GetElementalColor](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor.md)

# Called by

- [Apply](../../../../../functions/modeler/coloring/color_cells_in_touch/ColorCellsInTouch/Apply.md)