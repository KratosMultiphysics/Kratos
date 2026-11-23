---
type: C++ Method
title: ShapeFunctionsSecondDerivativesTransformOnAllIntegrationPoints
resource: utilities/geometry_utilities.cpp#L516-L538
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/GeometryUtils/ShapeFunctionsSecondDerivativesTransformOnIntegrationPoint
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/VerifyShapeFunctionsSecondDerivativesInterpolation
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GeometryUtils::ShapeFunctionsSecondDerivativesTransformOnAllIntegrationPoints( DenseVector<DenseVector<Matrix>>& rResult, const GeometryType& rGeometry, const GeometryType::IntegrationMethod& rIntegrationMethod )`

# Calls

- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [ShapeFunctionsSecondDerivativesTransformOnIntegrationPoint](../../../../functions/utilities/geometry_utilities/GeometryUtils/ShapeFunctionsSecondDerivativesTransformOnIntegrationPoint.md)

# Called by

- [VerifyShapeFunctionsSecondDerivativesInterpolation](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/VerifyShapeFunctionsSecondDerivativesInterpolation.md)