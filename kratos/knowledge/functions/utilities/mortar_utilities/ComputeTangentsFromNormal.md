---
type: C++ Function
title: ComputeTangentsFromNormal
resource: utilities/mortar_utilities.cpp#L210-L232
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/math_utils/OrthonormalBasis
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/mortar_utilities/ComputeNodesTangentFromNormalModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/mortar_utilities/ComputeTangentNodeWithLMAndSlip
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ComputeTangentsFromNormal( Node& rNode, const array_1d<double, 3>& rNormal, const std::size_t Dimension )`

# Calls

- [OrthonormalBasis](../../../functions/utilities/math_utils/OrthonormalBasis.md)

# Called by

- [ComputeNodesTangentFromNormalModelPart](../../../functions/utilities/mortar_utilities/ComputeNodesTangentFromNormalModelPart.md)
- [ComputeTangentNodeWithLMAndSlip](../../../functions/utilities/mortar_utilities/ComputeTangentNodeWithLMAndSlip.md)