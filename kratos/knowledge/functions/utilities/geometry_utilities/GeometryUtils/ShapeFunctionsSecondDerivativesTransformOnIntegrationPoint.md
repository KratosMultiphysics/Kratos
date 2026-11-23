---
type: C++ Method
title: ShapeFunctionsSecondDerivativesTransformOnIntegrationPoint
resource: utilities/geometry_utilities.cpp#L543-L689
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/geometry_utilities/GeometryUtils/ShapeFunctionsSecondDerivativesTransformOnAllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GeometryUtils::ShapeFunctionsSecondDerivativesTransformOnIntegrationPoint( const Matrix& DN_DX, const GeometryType& rGeometry, const GeometryType::CoordinatesArrayType& rLocalIntegrationPointCoordinates, DenseVector<Matrix>& rResult)`

# Calls

- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [prod](../../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [ShapeFunctionsSecondDerivativesTransformOnAllIntegrationPoints](../../../../functions/utilities/geometry_utilities/GeometryUtils/ShapeFunctionsSecondDerivativesTransformOnAllIntegrationPoints.md)