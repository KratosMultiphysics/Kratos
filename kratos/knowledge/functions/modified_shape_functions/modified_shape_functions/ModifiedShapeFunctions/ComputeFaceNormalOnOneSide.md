---
type: C++ Method
title: ComputeFaceNormalOnOneSide
resource: modified_shape_functions/modified_shape_functions.cpp#L688-L718
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveSideInterfaceAreaNormals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeNegativeSideInterfaceAreaNormals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveExteriorFaceAreaNormals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeNegativeExteriorFaceAreaNormals
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModifiedShapeFunctions::ComputeFaceNormalOnOneSide( AreaNormalsContainerType& rInterfaceAreaNormalValues, const std::vector<IndexedPointGeometryPointerType> &rInterfacesVector, const IntegrationMethodType IntegrationMethod)`

# Called by

- [ComputePositiveSideInterfaceAreaNormals](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveSideInterfaceAreaNormals.md)
- [ComputeNegativeSideInterfaceAreaNormals](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeNegativeSideInterfaceAreaNormals.md)
- [ComputePositiveExteriorFaceAreaNormals](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveExteriorFaceAreaNormals.md)
- [ComputeNegativeExteriorFaceAreaNormals](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeNegativeExteriorFaceAreaNormals.md)