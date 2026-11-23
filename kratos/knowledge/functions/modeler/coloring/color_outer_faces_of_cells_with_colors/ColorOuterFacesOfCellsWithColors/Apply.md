---
type: C++ Method
title: Apply
resource: modeler/coloring/color_outer_faces_of_cells_with_colors.cpp#L28-L51
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/coloring/color_outer_faces_of_cells_with_colors/ColorOuterFacesOfCellsWithColors/ApplyColorIfOuterFace
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ColorOuterFacesOfCellsWithColors::Apply() const`

# Calls

- [GetInt](../../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetElementalColor](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor.md)
- [ApplyColorIfOuterFace](../../../../../functions/modeler/coloring/color_outer_faces_of_cells_with_colors/ColorOuterFacesOfCellsWithColors/ApplyColorIfOuterFace.md)