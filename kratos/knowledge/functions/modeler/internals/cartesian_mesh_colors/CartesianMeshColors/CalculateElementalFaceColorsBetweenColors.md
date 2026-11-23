---
type: C++ Method
title: CalculateElementalFaceColorsBetweenColors
resource: modeler/internals/cartesian_mesh_colors.h#L555-L611
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/coloring/color_cell_faces_between_colors/ColorCellFacesBetweenColors/Apply
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateElementalFaceColorsBetweenColors(int FaceColor, int OutsideColor, int VolumeColor)`

# Calls

- [GetElementalColor](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor.md)

# Called by

- [Apply](../../../../../functions/modeler/coloring/color_cell_faces_between_colors/ColorCellFacesBetweenColors/Apply.md)