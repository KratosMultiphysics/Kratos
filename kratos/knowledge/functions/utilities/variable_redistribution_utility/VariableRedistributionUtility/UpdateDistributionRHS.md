---
type: C++ Method
title: UpdateDistributionRHS
resource: utilities/variable_redistribution_utility.cpp#L584-L619
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VariableRedistributionUtility::UpdateDistributionRHS( ModelPart& rModelPart, TContainerType& rEntitiesContainer, const Variable< TValueType >& rPointVariable, const Variable< TValueType >& rDistributedVariable, BoundedMatrix<double, TNumNodes, TNumNodes>& rMassMatrix)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)