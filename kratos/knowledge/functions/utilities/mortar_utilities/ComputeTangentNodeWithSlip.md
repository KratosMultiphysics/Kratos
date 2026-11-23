---
type: C++ Function
title: ComputeTangentNodeWithSlip
resource: utilities/mortar_utilities.cpp#L295-L357
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/OrthonormalBasis
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/mortar_utilities/ComputeNodesTangentModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ComputeTangentNodeWithSlip( Node& rNode, const std::size_t StepLM, const Variable<array_1d<double, 3>>* pSlipVariable, const double SlipCoefficient, const std::size_t Dimension )`

# Calls

- [norm_2](../../../functions/includes/amatrix_interface/norm_2.md)
- [OrthonormalBasis](../../../functions/utilities/math_utils/OrthonormalBasis.md)

# Called by

- [ComputeNodesTangentModelPart](../../../functions/utilities/mortar_utilities/ComputeNodesTangentModelPart.md)