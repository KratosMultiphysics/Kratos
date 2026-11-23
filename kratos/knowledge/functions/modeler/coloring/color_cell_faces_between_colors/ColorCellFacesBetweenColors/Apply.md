---
type: C++ Method
title: Apply
resource: modeler/coloring/color_cell_faces_between_colors.cpp#L28-L36
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateElementalFaceColorsBetweenColors
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ColorCellFacesBetweenColors::Apply() const`

# Calls

- [GetInt](../../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [CalculateElementalFaceColorsBetweenColors](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateElementalFaceColorsBetweenColors.md)