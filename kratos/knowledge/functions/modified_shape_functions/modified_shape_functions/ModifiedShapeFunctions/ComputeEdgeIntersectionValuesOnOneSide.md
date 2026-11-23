---
type: C++ Method
title: ComputeEdgeIntersectionValuesOnOneSide
resource: modified_shape_functions/modified_shape_functions.cpp#L721-L740
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeShapeFunctionsOnPositiveEdgeIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeShapeFunctionsOnNegativeEdgeIntersections
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModifiedShapeFunctions::ComputeEdgeIntersectionValuesOnOneSide( const Matrix &rPmatrix, Matrix &rEdgeShapeFunctionValues)`

# Calls

- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [ComputeShapeFunctionsOnPositiveEdgeIntersections](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeShapeFunctionsOnPositiveEdgeIntersections.md)
- [ComputeShapeFunctionsOnNegativeEdgeIntersections](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeShapeFunctionsOnNegativeEdgeIntersections.md)