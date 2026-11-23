---
type: C++ Method
title: VerifyShapeFunctionsSecondDerivativesValues
resource: utilities/geometry_tester.cpp#L964-L1025
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_frobenius
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestTetrahedra3D4N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestTriangle2D3N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestTriangle2D6N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestQuadrilateral2D4N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestQuadrilateral2D9N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestQuadrilateralInterface2D4N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestHexahedra3D8N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestHexahedra3D27N
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool GeometryTesterUtility::VerifyShapeFunctionsSecondDerivativesValues( GeometryType& rGeometry, GeometryType::CoordinatesArrayType& rGlobalCoordinates, std::stringstream& rErrorMessage ) const`

# Calls

- [norm_frobenius](../../../../functions/includes/amatrix_interface/norm_frobenius.md)

# Called by

- [StreamTestTetrahedra3D4N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestTetrahedra3D4N.md)
- [StreamTestTriangle2D3N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestTriangle2D3N.md)
- [StreamTestTriangle2D6N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestTriangle2D6N.md)
- [StreamTestQuadrilateral2D4N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestQuadrilateral2D4N.md)
- [StreamTestQuadrilateral2D9N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestQuadrilateral2D9N.md)
- [StreamTestQuadrilateralInterface2D4N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestQuadrilateralInterface2D4N.md)
- [StreamTestHexahedra3D8N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestHexahedra3D8N.md)
- [StreamTestHexahedra3D27N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestHexahedra3D27N.md)