---
type: C++ Method
title: ApplyColorToConnectedCellsInTouchWithGeometry
resource: modeler/coloring/color_connected_cells_in_touch.cpp#L59-L83
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
  - target: functions/modeler/coloring/color_connected_cells_in_touch/ColorConnectedCellsInTouch/ColorConnectedCellsToThisCell
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/coloring/color_connected_cells_in_touch/ColorConnectedCellsInTouch/Apply
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ColorConnectedCellsInTouch::ApplyColorToConnectedCellsInTouchWithGeometry( const Geometry<Node>& rGeometry, const int InsideColor, const int CellColor ) const`

# Calls

- [CalculateOuterMinMaxNodePositions](../../../../../functions/modeler/internals/cartesian_mesh_colors/CalculateOuterMinMaxNodePositions.md)
- [GetElementalColor](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor.md)
- [ColorConnectedCellsToThisCell](../../../../../functions/modeler/coloring/color_connected_cells_in_touch/ColorConnectedCellsInTouch/ColorConnectedCellsToThisCell.md)

# Called by

- [Apply](../../../../../functions/modeler/coloring/color_connected_cells_in_touch/ColorConnectedCellsInTouch/Apply.md)