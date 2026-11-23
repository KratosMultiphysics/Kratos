---
type: C++ Method
title: ComputeDomainSizeOnOneSide
resource: modified_shape_functions/modified_shape_functions.cpp#L742-L749
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveSideDomainSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeNegativeSideDomainSize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double ModifiedShapeFunctions::ComputeDomainSizeOnOneSide(const std::vector<IndexedPointGeometryPointerType> &rSubdivisionsVector) const`

# Called by

- [ComputePositiveSideDomainSize](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveSideDomainSize.md)
- [ComputeNegativeSideDomainSize](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeNegativeSideDomainSize.md)