---
type: C++ Method
title: GetModifiedShapeFunctionsValues
resource: utilities/embedded_skin_utility.cpp#L373-L398
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Matrix EmbeddedSkinUtility<TDim>::GetModifiedShapeFunctionsValues( const ModifiedShapeFunctions::UniquePointer &rpModifiedShapeFunctions, const std::string &rInterfaceSide) const`

# Calls

- [ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues.md)
- [ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues.md)