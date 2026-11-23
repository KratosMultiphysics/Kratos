---
type: C++ Method
title: GetNodalDistances
resource: utilities/divide_geometry.cpp#L106-L108
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modified_shape_functions/ausas_modified_shape_functions/AusasModifiedShapeFunctions/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/ausas_modified_shape_functions/AusasModifiedShapeFunctions/SetPositiveSideCondensationMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/ausas_modified_shape_functions/AusasModifiedShapeFunctions/SetNegativeSideCondensationMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/ausas_modified_shape_functions/AusasModifiedShapeFunctions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/tetrahedra_3d_4_ausas_incised_shape_functions/Tetrahedra3D4AusasIncisedShapeFunctions/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/tetrahedra_3d_4_ausas_incised_shape_functions/Tetrahedra3D4AusasIncisedShapeFunctions/SetPositiveSideCondensationMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/tetrahedra_3d_4_ausas_incised_shape_functions/Tetrahedra3D4AusasIncisedShapeFunctions/SetNegativeSideCondensationMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/tetrahedra_3d_4_ausas_incised_shape_functions/Tetrahedra3D4AusasIncisedShapeFunctions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions/Tetrahedra3D4AusasModifiedShapeFunctions/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions/Tetrahedra3D4AusasModifiedShapeFunctions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/tetrahedra_3d_4_modified_shape_functions/Tetrahedra3D4ModifiedShapeFunctions/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/tetrahedra_3d_4_modified_shape_functions/Tetrahedra3D4ModifiedShapeFunctions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/triangle_2d_3_ausas_incised_shape_functions/Triangle2D3AusasIncisedShapeFunctions/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/triangle_2d_3_ausas_incised_shape_functions/Triangle2D3AusasIncisedShapeFunctions/SetPositiveSideCondensationMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/triangle_2d_3_ausas_incised_shape_functions/Triangle2D3AusasIncisedShapeFunctions/SetNegativeSideCondensationMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/triangle_2d_3_ausas_incised_shape_functions/Triangle2D3AusasIncisedShapeFunctions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions/Triangle2D3AusasModifiedShapeFunctions/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions/Triangle2D3AusasModifiedShapeFunctions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/triangle_2d_3_modified_shape_functions/Triangle2D3ModifiedShapeFunctions/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/triangle_2d_3_modified_shape_functions/Triangle2D3ModifiedShapeFunctions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_tetrahedra_3d_4/DivideTetrahedra3D4<TPointType>/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_tetrahedra_3d_4/DivideTetrahedra3D4<TPointType>/GenerateDivision
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_tetrahedra_3d_4/DivideTetrahedra3D4
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_triangle_2d_3/DivideTriangle2D3<TPointType>/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_triangle_2d_3/DivideTriangle2D3<TPointType>/GenerateDivision
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_triangle_2d_3/DivideTriangle2D3
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_triangle_3d_3/DivideTriangle3D3
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Vector DivideGeometry<TPointType>::GetNodalDistances() const`

# Called by

- [PrintData](../../../../functions/modified_shape_functions/ausas_modified_shape_functions/AusasModifiedShapeFunctions/PrintData.md)
- [SetPositiveSideCondensationMatrix](../../../../functions/modified_shape_functions/ausas_modified_shape_functions/AusasModifiedShapeFunctions/SetPositiveSideCondensationMatrix.md)
- [SetNegativeSideCondensationMatrix](../../../../functions/modified_shape_functions/ausas_modified_shape_functions/AusasModifiedShapeFunctions/SetNegativeSideCondensationMatrix.md)
- [AusasModifiedShapeFunctions](../../../../functions/modified_shape_functions/ausas_modified_shape_functions/AusasModifiedShapeFunctions.md)
- [PrintData](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/PrintData.md)
- [PrintData](../../../../functions/modified_shape_functions/tetrahedra_3d_4_ausas_incised_shape_functions/Tetrahedra3D4AusasIncisedShapeFunctions/PrintData.md)
- [SetPositiveSideCondensationMatrix](../../../../functions/modified_shape_functions/tetrahedra_3d_4_ausas_incised_shape_functions/Tetrahedra3D4AusasIncisedShapeFunctions/SetPositiveSideCondensationMatrix.md)
- [SetNegativeSideCondensationMatrix](../../../../functions/modified_shape_functions/tetrahedra_3d_4_ausas_incised_shape_functions/Tetrahedra3D4AusasIncisedShapeFunctions/SetNegativeSideCondensationMatrix.md)
- [Tetrahedra3D4AusasIncisedShapeFunctions](../../../../functions/modified_shape_functions/tetrahedra_3d_4_ausas_incised_shape_functions/Tetrahedra3D4AusasIncisedShapeFunctions.md)
- [PrintData](../../../../functions/modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions/Tetrahedra3D4AusasModifiedShapeFunctions/PrintData.md)
- [Tetrahedra3D4AusasModifiedShapeFunctions](../../../../functions/modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions/Tetrahedra3D4AusasModifiedShapeFunctions.md)
- [PrintData](../../../../functions/modified_shape_functions/tetrahedra_3d_4_modified_shape_functions/Tetrahedra3D4ModifiedShapeFunctions/PrintData.md)
- [Tetrahedra3D4ModifiedShapeFunctions](../../../../functions/modified_shape_functions/tetrahedra_3d_4_modified_shape_functions/Tetrahedra3D4ModifiedShapeFunctions.md)
- [PrintData](../../../../functions/modified_shape_functions/triangle_2d_3_ausas_incised_shape_functions/Triangle2D3AusasIncisedShapeFunctions/PrintData.md)
- [SetPositiveSideCondensationMatrix](../../../../functions/modified_shape_functions/triangle_2d_3_ausas_incised_shape_functions/Triangle2D3AusasIncisedShapeFunctions/SetPositiveSideCondensationMatrix.md)
- [SetNegativeSideCondensationMatrix](../../../../functions/modified_shape_functions/triangle_2d_3_ausas_incised_shape_functions/Triangle2D3AusasIncisedShapeFunctions/SetNegativeSideCondensationMatrix.md)
- [Triangle2D3AusasIncisedShapeFunctions](../../../../functions/modified_shape_functions/triangle_2d_3_ausas_incised_shape_functions/Triangle2D3AusasIncisedShapeFunctions.md)
- [PrintData](../../../../functions/modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions/Triangle2D3AusasModifiedShapeFunctions/PrintData.md)
- [Triangle2D3AusasModifiedShapeFunctions](../../../../functions/modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions/Triangle2D3AusasModifiedShapeFunctions.md)
- [PrintData](../../../../functions/modified_shape_functions/triangle_2d_3_modified_shape_functions/Triangle2D3ModifiedShapeFunctions/PrintData.md)
- [Triangle2D3ModifiedShapeFunctions](../../../../functions/modified_shape_functions/triangle_2d_3_modified_shape_functions/Triangle2D3ModifiedShapeFunctions.md)
- [PrintData](../../../../functions/utilities/divide_tetrahedra_3d_4/DivideTetrahedra3D4<TPointType>/PrintData.md)
- [GenerateDivision](../../../../functions/utilities/divide_tetrahedra_3d_4/DivideTetrahedra3D4<TPointType>/GenerateDivision.md)
- [DivideTetrahedra3D4](../../../../functions/utilities/divide_tetrahedra_3d_4/DivideTetrahedra3D4.md)
- [PrintData](../../../../functions/utilities/divide_triangle_2d_3/DivideTriangle2D3<TPointType>/PrintData.md)
- [GenerateDivision](../../../../functions/utilities/divide_triangle_2d_3/DivideTriangle2D3<TPointType>/GenerateDivision.md)
- [DivideTriangle2D3](../../../../functions/utilities/divide_triangle_2d_3/DivideTriangle2D3.md)
- [DivideTriangle3D3](../../../../functions/utilities/divide_triangle_3d_3/DivideTriangle3D3.md)