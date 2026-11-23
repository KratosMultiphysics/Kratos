---
type: C++ Method
title: ApplyColorIfOuterFace
resource: modeler/coloring/color_outer_faces_of_cells_with_colors.cpp#L54-L84
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/coloring/color_outer_faces_of_cells_with_colors/ColorOuterFacesOfCellsWithColors/Apply
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ColorOuterFacesOfCellsWithColors::ApplyColorIfOuterFace( const int InterfaceColor, const int CellColor, const array_1d<std::size_t, 3>& rCellIndices) const`

# Calls

- [GetElementalColor](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor.md)

# Called by

- [Apply](../../../../../functions/modeler/coloring/color_outer_faces_of_cells_with_colors/ColorOuterFacesOfCellsWithColors/Apply.md)