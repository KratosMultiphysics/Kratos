---
type: C++ Method
title: ColorConnectedCellsToThisCell
resource: modeler/coloring/color_connected_cells_in_touch.cpp#L86-L137
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/pointer_vector/emplace_back
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/pointer_vector_set/PointerVectorSet/pop_back
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/coloring/color_connected_cells_in_touch/ColorConnectedCellsInTouch/ApplyColorToConnectedCellsInTouchWithGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ColorConnectedCellsInTouch::ColorConnectedCellsToThisCell( const std::size_t I, const std::size_t J, const std::size_t K, const int InsideColor, const int CellColor, CartesianMeshColors& rColors) const`

# Calls

- [emplace_back](../../../../../functions/containers/pointer_vector/emplace_back.md)
- [pop_back](../../../../../functions/containers/pointer_vector_set/PointerVectorSet/pop_back.md)
- [GetElementalColor](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor.md)

# Called by

- [ApplyColorToConnectedCellsInTouchWithGeometry](../../../../../functions/modeler/coloring/color_connected_cells_in_touch/ColorConnectedCellsInTouch/ApplyColorToConnectedCellsInTouchWithGeometry.md)