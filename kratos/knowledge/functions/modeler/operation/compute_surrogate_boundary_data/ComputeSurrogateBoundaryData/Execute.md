---
type: C++ Method
title: Execute
resource: modeler/operation/compute_surrogate_boundary_data.cpp#L40-L133
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/InitializeRays
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateNodalRayColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/operation/compute_surrogate_boundary_data/SurrogateBoundaryNode/SetNodePointer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/spatial_search_result/SpatialSearchResult/IsObjectFound
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/operation/compute_surrogate_boundary_data/ComputeSurrogateBoundaryData/PrintSurrogateBoundaryData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ComputeSurrogateBoundaryData::Execute()`

# Calls

- [GetInt](../../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetString](../../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [InitializeRays](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/InitializeRays.md)
- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [CalculateNodalRayColors](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateNodalRayColors.md)
- [SetNodePointer](../../../../../functions/modeler/operation/compute_surrogate_boundary_data/SurrogateBoundaryNode/SetNodePointer.md)
- [IsObjectFound](../../../../../functions/spatial_containers/spatial_search_result/SpatialSearchResult/IsObjectFound.md)
- [norm_2](../../../../../functions/includes/amatrix_interface/norm_2.md)
- [PrintSurrogateBoundaryData](../../../../../functions/modeler/operation/compute_surrogate_boundary_data/ComputeSurrogateBoundaryData/PrintSurrogateBoundaryData.md)