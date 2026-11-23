---
type: C++ Method
title: Apply
resource: modeler/coloring/color_cells_with_center_in_sphere_around_nodes.cpp#L29-L59
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetCenterOfElement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ColorCellsWithCenterInSphereArounNodes::Apply() const`

# Calls

- [GetInt](../../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetCenterOfElement](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetCenterOfElement.md)
- [GetElementalColor](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor.md)