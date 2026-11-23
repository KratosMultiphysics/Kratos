---
type: C++ Method
title: ComputeShapeFunctionsOnNegativeEdgeIntersections
resource: modified_shape_functions/modified_shape_functions.cpp#L427-L441
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeEdgeIntersectionValuesOnOneSide
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/GetModifiedShapeFunctionsValuesOnEdge
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModifiedShapeFunctions::ComputeShapeFunctionsOnNegativeEdgeIntersections(Matrix &rNegativeEdgeIntersectionsShapeFunctionsValues)`

# Calls

- [ComputeEdgeIntersectionValuesOnOneSide](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeEdgeIntersectionValuesOnOneSide.md)

# Called by

- [GetModifiedShapeFunctionsValuesOnEdge](../../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/GetModifiedShapeFunctionsValuesOnEdge.md)