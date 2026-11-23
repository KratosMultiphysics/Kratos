---
type: C++ Function
title: CalculateMinMaxCenterOfElementPositions
resource: modeler/internals/cartesian_mesh_colors.h#L663-L681
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/coloring/color_cell_faces/ColorCellFaces/Apply
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/coloring/color_cells_with_inside_center/ColorCellsWithInsideCenter/Apply
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/AddGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateMinMaxCenterOfElementPositions(TPointsContainerType const& Points, array_1d< std::size_t, 3 >& MinNodePosition, array_1d< std::size_t, 3 >& MaxNodePosition)`

# Called by

- [Apply](../../../../functions/modeler/coloring/color_cell_faces/ColorCellFaces/Apply.md)
- [Apply](../../../../functions/modeler/coloring/color_cells_with_inside_center/ColorCellsWithInsideCenter/Apply.md)
- [AddGeometry](../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/AddGeometry.md)