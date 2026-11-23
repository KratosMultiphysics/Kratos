---
type: C++ Method
title: ComputeValuesOnOneSide
resource: modified_shape_functions/modified_shape_functions.cpp#L487-L571
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModifiedShapeFunctions::ComputeValuesOnOneSide( Matrix &rShapeFunctionsValues, ShapeFunctionsGradientsType &rShapeFunctionsGradientsValues, Vector &rWeightsValues, const std::vector<IndexedPointGeometryPointerType> &rSubdivisionsVector, const Matrix &rPmatrix, const IntegrationMethodType IntegrationMethod)`

# Calls

- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [Weight](../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)
- [prod](../../../../functions/includes/amatrix_interface/prod.md)