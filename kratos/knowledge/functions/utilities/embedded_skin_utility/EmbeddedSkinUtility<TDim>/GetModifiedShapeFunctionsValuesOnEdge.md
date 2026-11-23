---
type: C++ Method
title: GetModifiedShapeFunctionsValuesOnEdge
resource: utilities/embedded_skin_utility.cpp#L401-L415
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeShapeFunctionsOnPositiveEdgeIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeShapeFunctionsOnNegativeEdgeIntersections
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/embedded_skin_utility/InterpolateMeshVariableToSkinSpecialization
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Matrix EmbeddedSkinUtility<TDim>::GetModifiedShapeFunctionsValuesOnEdge( const ModifiedShapeFunctions::UniquePointer &rpModifiedShapeFunctions, const std::string &rInterfaceSide) const`

# Calls

- [ComputeShapeFunctionsOnPositiveEdgeIntersections](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeShapeFunctionsOnPositiveEdgeIntersections.md)
- [ComputeShapeFunctionsOnNegativeEdgeIntersections](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeShapeFunctionsOnNegativeEdgeIntersections.md)

# Called by

- [InterpolateMeshVariableToSkinSpecialization](../../../../functions/utilities/embedded_skin_utility/InterpolateMeshVariableToSkinSpecialization.md)