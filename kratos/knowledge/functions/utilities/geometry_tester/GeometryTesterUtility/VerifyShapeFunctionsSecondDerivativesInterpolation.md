---
type: C++ Method
title: VerifyShapeFunctionsSecondDerivativesInterpolation
resource: utilities/geometry_tester.cpp#L1027-L1105
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/GeometryUtils/ShapeFunctionsSecondDerivativesTransformOnAllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestQuadrilateral2D9N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestHexahedra3D27N
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool GeometryTesterUtility::VerifyShapeFunctionsSecondDerivativesInterpolation( GeometryType& rGeometry, std::stringstream& rErrorMessage ) const`

# Calls

- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [ShapeFunctionsSecondDerivativesTransformOnAllIntegrationPoints](../../../../functions/utilities/geometry_utilities/GeometryUtils/ShapeFunctionsSecondDerivativesTransformOnAllIntegrationPoints.md)

# Called by

- [StreamTestQuadrilateral2D9N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestQuadrilateral2D9N.md)
- [StreamTestHexahedra3D27N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestHexahedra3D27N.md)