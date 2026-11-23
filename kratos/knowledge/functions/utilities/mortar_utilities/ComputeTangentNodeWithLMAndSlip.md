---
type: C++ Function
title: ComputeTangentNodeWithLMAndSlip
resource: utilities/mortar_utilities.cpp#L237-L290
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/mortar_utilities/ComputeTangentsFromNormal
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/mortar_utilities/ComputeNodesTangentModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ComputeTangentNodeWithLMAndSlip( Node& rNode, const std::size_t StepLM, const Variable<array_1d<double, 3>>* pSlipVariable, const double SlipCoefficient, const std::size_t Dimension )`

# Calls

- [IsNot](../../../functions/containers/flags/IsNot.md)
- [norm_2](../../../functions/includes/amatrix_interface/norm_2.md)
- [ComputeTangentsFromNormal](../../../functions/utilities/mortar_utilities/ComputeTangentsFromNormal.md)

# Called by

- [ComputeNodesTangentModelPart](../../../functions/utilities/mortar_utilities/ComputeNodesTangentModelPart.md)