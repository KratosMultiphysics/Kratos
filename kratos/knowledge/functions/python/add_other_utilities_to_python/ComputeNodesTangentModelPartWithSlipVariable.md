---
type: C++ Function
title: ComputeNodesTangentModelPartWithSlipVariable
resource: python/add_other_utilities_to_python.cpp#L114-L122
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/mortar_utilities/ComputeNodesTangentModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ComputeNodesTangentModelPartWithSlipVariable( ModelPart& rModelPart, const Variable<array_1d<double, 3>>& rSlipVariable, const double SlipCoefficient, const bool SlipAlways )`

# Calls

- [ComputeNodesTangentModelPart](../../../functions/utilities/mortar_utilities/ComputeNodesTangentModelPart.md)