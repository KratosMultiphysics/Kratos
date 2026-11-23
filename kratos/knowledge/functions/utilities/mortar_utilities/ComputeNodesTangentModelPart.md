---
type: C++ Function
title: ComputeNodesTangentModelPart
resource: utilities/mortar_utilities.cpp#L149-L180
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/mortar_utilities/ComputeTangentNodeWithLMAndSlip
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/mortar_utilities/ComputeTangentNodeWithSlip
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/add_other_utilities_to_python/ComputeNodesTangentModelPartWithSlipVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_other_utilities_to_python/ComputeNodesTangentModelPartWithOutSlipVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_other_utilities_to_python/ComputeNodesTangentModelPartWithSlipVariableNotAlwaysSlip
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_other_utilities_to_python/ComputeNodesTangentModelPartWithOutSlipVariableNotAlwaysSlip
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_other_utilities_to_python/ComputeNodesTangentModelPartWithSlipVariableNotAlwaysSlipUnitary
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_other_utilities_to_python/ComputeNodesTangentModelPartWithOutSlipVariableNotAlwaysSlipUnitary
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ComputeNodesTangentModelPart( ModelPart& rModelPart, const Variable<array_1d<double, 3>>* pSlipVariable, const double SlipCoefficient, const bool SlipAlways )`

# Calls

- [ComputeTangentNodeWithLMAndSlip](../../../functions/utilities/mortar_utilities/ComputeTangentNodeWithLMAndSlip.md)
- [ComputeTangentNodeWithSlip](../../../functions/utilities/mortar_utilities/ComputeTangentNodeWithSlip.md)

# Called by

- [ComputeNodesTangentModelPartWithSlipVariable](../../../functions/python/add_other_utilities_to_python/ComputeNodesTangentModelPartWithSlipVariable.md)
- [ComputeNodesTangentModelPartWithOutSlipVariable](../../../functions/python/add_other_utilities_to_python/ComputeNodesTangentModelPartWithOutSlipVariable.md)
- [ComputeNodesTangentModelPartWithSlipVariableNotAlwaysSlip](../../../functions/python/add_other_utilities_to_python/ComputeNodesTangentModelPartWithSlipVariableNotAlwaysSlip.md)
- [ComputeNodesTangentModelPartWithOutSlipVariableNotAlwaysSlip](../../../functions/python/add_other_utilities_to_python/ComputeNodesTangentModelPartWithOutSlipVariableNotAlwaysSlip.md)
- [ComputeNodesTangentModelPartWithSlipVariableNotAlwaysSlipUnitary](../../../functions/python/add_other_utilities_to_python/ComputeNodesTangentModelPartWithSlipVariableNotAlwaysSlipUnitary.md)
- [ComputeNodesTangentModelPartWithOutSlipVariableNotAlwaysSlipUnitary](../../../functions/python/add_other_utilities_to_python/ComputeNodesTangentModelPartWithOutSlipVariableNotAlwaysSlipUnitary.md)