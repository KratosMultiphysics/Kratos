---
type: C++ Method
title: Apply
resource: modeler/coloring/color_cells_with_inside_center.cpp#L28-L66
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CalculateMinMaxCenterOfElementPositions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/InitializeRays
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateElementalRayColors
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ColorCellsWithInsideCenter::Apply() const`

# Calls

- [GetString](../../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetInt](../../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [CalculateMinMaxCenterOfElementPositions](../../../../../functions/modeler/internals/cartesian_mesh_colors/CalculateMinMaxCenterOfElementPositions.md)
- [InitializeRays](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/InitializeRays.md)
- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [CalculateElementalRayColors](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateElementalRayColors.md)