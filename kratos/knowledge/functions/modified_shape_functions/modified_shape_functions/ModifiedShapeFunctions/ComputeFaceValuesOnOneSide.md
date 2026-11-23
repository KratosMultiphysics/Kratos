---
type: C++ Method
title: ComputeFaceValuesOnOneSide
resource: modified_shape_functions/modified_shape_functions.cpp#L574-L685
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
  called_by:
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModifiedShapeFunctions::ComputeFaceValuesOnOneSide( Matrix &rInterfaceShapeFunctionsValues, ShapeFunctionsGradientsType &rInterfaceShapeFunctionsGradientsValues, Vector &rInterfaceWeightsValues, const std::vector<IndexedPointGeometryPointerType> &rInterfacesVector, const std::vector<IndexedPointGeometryPointerType> &rParentGeometriesVector, const std::vector<unsigned int> &rInterfacesParentIdsVector, const Matrix &rPmatrix, const IntegrationMethodType IntegrationMethod)`

# Calls

- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [Weight](../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)
- [prod](../../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues.md)
- [ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues.md)
- [ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues.md)
- [ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues.md)